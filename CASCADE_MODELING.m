%%%% CASCADE Master script. This script contains all necessary sub-scripts
%%%% and functions to either 
% a) Run the model from the raw data
% b) Load data from my previous experiments and look only at specific parts
% of the modeling process
% c) Make runs of the CASCADE model and look at its results. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define output folder for the preprocessed data
PreprocessedOutFolder='CASCADE_PreProcessing_Output'; 

% This was added as a test for GitHub

addpath(genpath(pwd))

% Run the master pre-processing from raw data, or load preprocessed data from a file 
    prompt = 'Do you want to rerun preprocessing? [Y/N]: ';
    str = input(prompt,'s');

    if strcmpi('y',str)
        
       % load raw date 
       Import_data_20150911
       global outlet_node s; outlet_node=48865;  
       
       %run preprocessing
       CASCADE_PREPROCESSING
       
    else 
%         filetoload='CASCADE_PreProcessingOutput_20160111_1500';
          filetoload='201509122_temp';
            prompt = ['Do you want to load the file:  ' filetoload '  [Y/N]: '];
            str = input(prompt,'s');           
            
            if strcmpi('y',str)
                load(filetoload) 
                                            
                if exist('Network')~=1
                    MultiGraphDefinition
                end
                
            end 
           
    end; clear input
 
    global outlet_node_new
%%%%%%%% Make some figures for quality control %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%SlpClasses=[1e-5 1e-4 1e-3  1e-2  1e-1 1]
%networkPlotter(AggData,AggData(:,ID_Slp),SlpClasses, AggData(:,ID_StrO),'parula','Slope','Slope  [m/m]',0);
%networkPlotter(raw_data_sort,dissolve,SlpClasses, raw_data_sort(:,ID_StrO),'parula','Reach Aggregation',[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate hydraulic parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('hydraulicData','var')==1
prompt = 'Do you want to run hydraulic calculations? [Y/N]: ';
str = input(prompt,'s');
elseif  exist('hydraulicData','var')==0 || strcmpi('y',str)==1
hydraulicData=hydraulicCalc(AggData, 1); % ATTENTION: hydraulicCalc was slightly changed in comparison to the previous version. 
    % Hydraulic data contains: 1: Water level; 2: local grainsize; 3: kstrickler; 4: flow velocity; 5: Froude; 6: critical shear stress ; % Output storage for hydrological variables 

    
    
        %%%%%%%% Make some figures for quality control %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        hClasses=(2:5:80); % Water level at bankfull
        networkPlotter(AggData,hydraulicData(:,1),hClasses, AggData(:,ID_StrO),'parula','Water level','Water level  [m]',0);

        FRClasses=(0.1:0.2:3); % Froude at bankfull
        networkPlotter(AggData,hydraulicData(:,5),FRClasses, AggData(:,ID_StrO),'parula','Froude','FR',0);

        vClasses=(1:5:100); % velocity at bankfull
        networkPlotter(AggData,hydraulicData(:,4),vClasses, AggData(:,ID_StrO),'parula','Velocity','v [m/s]',0);

        kstClasses=(10:5:50); % velocity at bankfull
        networkPlotter(AggData,hydraulicData(:,3),kstClasses, AggData(:,ID_StrO),'parula','k strickler','kst [m/s^1/3]',0);

end

        % Key plots: spatial and probability distribution of grain size.
        figure('name','Grain size distribution')
            subplot(1,2,1)
                scatter(AggData(:,ID_FX), AggData(:,ID_FY),(hydraulicData(:,2)*50),((hydraulicData(:,2))*1000),'filled')
                cb=colorbar;  

            subplot(1,2,2)
                Wentworth=[0.125 0.25 0.5 1 2 4 8 16 32 64 256 512]*1E-3;
                bar(histcounts(hydraulicData(:,2),Wentworth))
                set(gca,'XTickLabels',(num2str(Wentworth')));
                xlabel('Grainsize [mm]'),ylabel('f [ ]')
                
 Dmat=parameterToMatrix(hydraulicData(:,2),II,2); % write Dmat              
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate QS_e^s for all cascades 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  prompt = 'Do you want to re-calculate the reference tranport capacity ? [Y/N]: ';
  str = input(prompt,'s');

  
  if strcmpi('y',str) % recalculate reference transport capacity 
% prepare some matrix inputs 
    dii=hydraulicData(:,2);
    Dmat=parameterToMatrix(dii,II,2);
    Wmat=parameterToMatrix(AggData(:,ID_Wac),II,1);
    Lmat=parameterToMatrix(AggData(:,ID_Length),II,1);

    HydraulicCalcs % -> result is TranspC_mat1, containing the NOT competition corrected transport capacity       

  else % load reference transport capacity from file
      
      filetoload='Transport_Capacity_Outputs_20160113_1304';
      prompt = ['Do you want to load ' filetoload ' the reference tranport capacity ? [Y/N]: '];
      str = input(prompt,'s');
      if strcmpi('y',str)
          load(filetoload)
          disp(['Loading: ' filetoload])
      end
      
  end
  
  TranspC_mat1=Qsij; % Define the original transport capacity matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make a single run of CASCADE to define connectivity in an undisturbed state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define some options: 
          scenario='Scenario 3'; 
          Competition=3; % 1: no competition, transport capacity is infinite (all inputs transported downstream), 2: no competition, but local transport capacity limits, 3: competition between pathways, and local capacity limits 
          Tcap_redistribution=0; % 1.: free trasnport capacity can be redistributed between pathways 
          transport_treshold=1e-2; % after how much percent of the initial input a pathway is interrupted 
          maxSupplyIncrease = 0.0; % select if and how much more sources downstream of reservoirs can produce
          addReservoir='No'; % use the following lines to add reservoirs
          reservoirFromN=[2361  2661  2881];
    % Run the CASCADE 
           calculationOrderType='dod';
           CASCADE_v6_fast
           Output(Output==0)=nan;
           disp(['TOTAL BASIN SEDIMET OUTPUT IS ' num2str(nansum(Input(:,outlet_node_new(1))),'%10.2e\n')])
           
        %%%% store some outputs as struct for easyly passing them to
        %%%% functions etc.
        
        global CascadeResults

         CascadeResults.Input=Input; 
         CascadeResults.Output=Output; 
         CascadeResults.Delivery=PP4; 
         CascadeResults.Network=Network; 
         CascadeResults.TransportCapacity=TranspC_mat2; 
         CascadeResults.Dmat=Dmat;
        
%%%%% Make some additional plots to check results

        %%%% Disconnectivity plot
        networkPlotter(AggData,nan(length(AggData),1),1, AggData(:,ID_StrO),'parula','Disconnectivity bottlenecks',[],0);
        hold on
        PlotDisconnectivity
    
        %%% Sediment transport from node 1 to the outlet
%         PLOT_sediment_transport_node_1_to_out_v2
    
        %%% Sediment delivery to the outlet 
        valOutput=find(~isnan(Output(:,outlet_node_new(1))) & Output(:,outlet_node_new(1))>0); 
        deliveryClasses=prctile(Output(valOutput,outlet_node_new(1)),5:5:99);
        deliveryClasses=[1e6 5e6 1e7 5e7 1e8 5e8 1e9 5e9 1e10];
        h_net=networkPlotter(AggData,Output(:,outlet_node_new(1)),deliveryClasses, AggData(:,ID_StrO),'parula','Delivery to basin outlet','Transport to outlet [kg/yr]',0);

        
        %%% Travel distance
                            for iii=1:length(Output); 
                               distToOut(iii)=Network.Downstream.Distance{iii}(outlet_node_new(1)); % total distance from outlet
                               Path_to_out=Network.Downstream.Path{iii}{outlet_node_new(1)}; % find the path to the outlet
                               if isempty(Path_to_out) Path_to_out=Network.Downstream.Path{iii}{outlet_node_new(2)}; end 
                               interruptionNode=find(isnan(CascadeResults.Output(iii,Path_to_out)),1,'first'); % find the first node where the path is interrupted
                               if isempty(interruptionNode) interruptionNode=length(Path_to_out); end % if the pathway is not interrupted plot tranport to the basin outlet    
                               interruptionNodeID=Path_to_out(interruptionNode); % absolute ID of interruption node
                               distToInterruption(iii)=Network.Downstream.Distance{iii}(interruptionNodeID)./1000;
                            end    
        distanceClass=prctile(distToInterruption,5:5:99);distanceClass(distanceClass<=0)=[];
        h_net=networkPlotter(AggData, distToInterruption,distanceClass, AggData(:,ID_StrO),'parula','INTERACTIVE PLOT. CLICK ON THE RIVER NETWORK','Delivery distance [km]',0,'interactive');
