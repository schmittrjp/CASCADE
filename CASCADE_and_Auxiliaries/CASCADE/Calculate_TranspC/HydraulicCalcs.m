%% Calculates the potential trnaport capacity (hence if there woould be no other grains transported) in all downstream nodes. 
% requires hydraulic and hydologic data as input 
% 3 Step procedure: 
% 1: Load observed hydrographs 
% 2: Hydrologic and hydraulic analyses: Derive a downscaled hydrograph for
% each reach. Dissect hydrpgraph in percentiles. Calculate hydraulic
% conditions in each percentile. 
% 3: Calculate QSij in each percentile 

%%  1: Load observed hydrographs 
%%% load discharge observations 

    directory_name = 'E:\Da River GIS Workspace\Hydrologic modelling\Matlab spatial calculations\Hydrologic observations';
    files = dir(directory_name);
    Q_data= cell(4,length(files)-2); % prepare storage array (first two entries in the files are void). 
    for ff=4:size(files,1)
        load([directory_name '\' files(ff,1).name])
        
    %%% calculate normalized discharge values, noamilzed by drainage area 
    station_name{ff}=files(ff,1).name;  station_name{ff}= station_name{ff}(1:end-4) % get the station name (without '.mat')
    
    Q_data{1,ff}=station_name{ff}; % get sttaion name 
    Q_data{2,ff}=getfield(eval(station_name{ff}),'StationID'); % get staion ID 
    % get Q data
    Q_data{3,ff}=eval([station_name{ff} '.signal.value']);
    Q_data{3,ff}=eval([station_name{ff} '.signal.value']); 
    %calculate Q1.5 for each station 
    pQ_1_5=(1-1/(1.5*365))*100; % non exceedance probability of Q1.5
    Q_data{4,ff}=prctile(Q_data{3,ff}, pQ_1_5);
%   Q_data{5,ff}=1.321*.^0.8214
    end

    Q_data(:,1:3)=[];
   
%% 2: Hydrologic and hydraulic analyses. Get the flow percentile values for all reaches and calculate the relevant hydraulic parameters at these percentile values.  All outputs are store in the statistics  array "stat". 
%%% Row-wise, stat contains 
%1.: Hydrograph scaling factor J,  J describes by how much to downscale the local hydrograph.
%2.: Percentile values either in fixed steps or according to a normal distribution 
%3.: The mean Q within each percentile "bin"
%4.: The number of observations in each percentile "bin"
%5.: Mean waterlevel in each percentile "bin"
%6.: Mean velocity in each percentile bin 
%7.: Total number of days on record  

% define DEM cell size 
   cellsize=30*30; % dem cellsize
% Define percentile values to be used 
   prc_val=0:5:100; % equally spaced (all percentiles conatin the same number of values)
   prc_val=[0.1 2.3 15.9 50 84.1 97.7 99.9]; % spaced according to 1 stdev in standard distribution (http://upload.wikimedia.org/wikipedia/en/5/5c/PR_and_NCE.gif)
% initialize stats 
    stats=cell(5,size(AggData,1));

% define options for hydraulic solving procedure

global wac q15 slp d50  kst rho s fminout kk taucrit
options = optimset;options.Algorithm = 'interior-point';options.GradObj = 'off';options.Hessian = 'lbfgs';options.Display = 'off';options.MaxIter = 50000;options.MaxFunEvals = 50000;options.TolFun = 0.00000000001;options.TolX = 0.00000000001;option.OutputFcn=@outfun;

% % load the complete, corrected Wac data
% load('Complete_Wac.mat')
% [~, Wac_resort_ind]=sort(raw_data(:,ID_FromN)); % WAC is not sorted, get the order of From nodes from the raw data
% corr_Wac=corr_Wac(Wac_resort_ind); clear Wac_resort_ind % resort WAc
  
multiWaitbar('Calculating percentile hydraulic conditions' ,0)
% multiWaitbar('Calculating percentile sediment transport' ,0)


    for ii=1:size(AggData,1) % loop through all reaches 
    multiWaitbar('Calculating percentile hydraulic conditions' ,ii/length(AggData)) ;      
        if AggData(ii,ID_SubWS)~= 0 % check if the current reach is within any ctachment 
        q15=(1.321*(AggData(ii,ID_Ad)).^0.8214); % local Q1.5, calculated with formula from MSc thesis, correct Ad to be in km2
        cid=AggData(ii,ID_SubWS) ; % catchment id of current reach
            cpos=find(cell2mat(Q_data(2,:))==cid); % where is this station stored in Q_data?                
            
    %%% calculate and store J 
            j_temp= q15/Q_data{4,cpos}; % Calculate the hydrograph scaling factor
            stats{1,AggData(ii,ID_FromN)}= j_temp; % Save J, ordered by the the Fromnode         
            
            % calculate and store the prctile values
            Q_reach=Q_data{3,cpos}*stats{1,AggData(ii,ID_FromN)}; %estimate the local hydrograph by multiplying the observation with teh scale factor
            stats{2,AggData(ii,ID_FromN)}= prctile(Q_reach,prc_val); % save the precentile values 
            stats{3,AggData(ii,ID_FromN)}=(stats{2,AggData(ii,ID_FromN)}(2:end)+stats{2,AggData(ii,ID_FromN)}(1:end-1))./2;   % average Q within a percentile  
           
            
            clear obs_prc
            for prc=2:length(stats{2,AggData(ii,ID_FromN)}) % loop through the percentile values of the current reach
                
                stats{4,AggData(ii,ID_FromN)}(prc)=sum(Q_reach>=stats{2,AggData(ii,ID_FromN)}(prc-1) & Q_reach<stats{2,AggData(ii,ID_FromN)}(prc)); % find how many values fall between two percentiles
            end
            
 
   %%% Calculate the hydraulic condictions for EACH percentile in EACH reach        
           
           %%% get the parameter that are constant for a reach 
           wac=AggData(ii,ID_Wac);
           slp=AggData(ii,ID_Slp);
           d50=mean(Dmat(Dmat(:,AggData(ii,ID_FromN))>0,AggData(ii,ID_FromN)));
            %define initial conditions
           h_init=wac/1;
            
       %%% Loop through all Q percentile values for that reach
            for prc=1:length(stats{3,AggData(ii,ID_FromN)}) 
           
            q15=stats{3,AggData(ii,ID_FromN)}(prc); % get the mean Q in the current bin 
            
%                 if prc==length(stats{3,AggData(ii,ID_FromN)}) % check if the flow in the largest percentile is smaller than the Q1.5. This can happen because of teh differnet methods used to derive teh Q1.5 (Ad scaling) versus the derivation of from the downscaled hydrograph
%                     if q15<raw_data_sort(raw_data2(ii,4),18) 
%                         q15=raw_data_sort(raw_data2(ii,4),18); % if this is the case, replace the flow in the largest percentile by the Q1.5 form the Ad interpolation. 
%                     end 
%                 end 

            %define initial conditions
            h_init=wac/2;
            %solve for h
            kk=1;
            [h,fval,exitflag,output]=fmincon(@hydraulic_solver_known_d50,h_init,[],[],[],[],0.2,20,[],options);
            v=q15/(h*wac);
            if v>20
            asd=212;
            end 
            
            stats{5,AggData(ii,ID_FromN)}(prc)=h;
            stats{6,AggData(ii,ID_FromN)}(prc)=v;
            end  
            
        end
stats{7,AggData(ii,ID_FromN)}=length(Q_data{3,cpos}); % total number of days on the record

    end
 
%% calculate reference transport capacity (i.e., TranspC_mat1)    
 CalculateQsij_par

%% save outputs to a file

    prompt = 'Do you want to store results of Transport capacity calculations [Y/N]: '
    str = input(prompt,'s');

if strcmpi('Y',str) 
       
    % check if top-level output folder exists 
       if exist(PreprocessedOutFolder,'dir')~=7 % if it doesn't exist, make and output folder
        mkdir(PreprocessedOutFolder)    
       end
    % check if daily output folder already exists 
       if exist(['Output_' datestr(now,'yyyymmdd')],'dir')~=7 % if it doesn't exist, make and output folder
         folder_name=['Output_' datestr(now,'yyyymmdd')]; % define subfolder name 
         mkdir(PreprocessedOutFolder,['Output_' datestr(now,'yyyymmdd')]) % make the daily output folder 
         mkdir([PreprocessedOutFolder '/' ['Output_' datestr(now,'yyyymmdd')]],'Matfiles') % make the output folder for the matfiles
         tCapOutputFolder=[PreprocessedOutFolder '/' ['Output_' datestr(now,'yyyymmdd')] '/' 'Matfiles']; % create a variable containing the name of the new output folder, for later use
         addpath(genpath([PreprocessedOutFolder '/' ['Output_' datestr(now,'yyyymmdd')]])) % add the folder to the search path
       end

    %store m-file

        save([tCapOutputFolder '/' 'Transport_Capacity_Outputs_' datestr(now,'yyyymmdd_HHMM')],'stats','QSij','Intermittency_mat') 

    
end

% 
% QS_classes=([ 1E4 1E5 1E6 1E7 1E8 1E9 2);
% networkPlotter(AggData,meanQSij,QS_classes, AggData(:,ID_StrO),'parula','Q_S Classes','Q_S [mio t / yr ]',0)

    
