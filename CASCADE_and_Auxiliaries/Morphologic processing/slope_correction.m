function [ raw_data_sort,SS_cor_results ] = slope_correction( raw_data_sort,SSCorMeth,D )
% Corrects the slopes measured from a DEM. Specifically, it repairs slope
% for river reaches in which slope information was lost during filling
% operations.

% There are three slope correction methods methods

% Method 1 ('SSCorMeth 1'), adapts the slope of the void reach, and the
% slope of the downstream reaches, in order to derive a new slope for both.
%->SLOPE IS A Function of measured morphologic conditions, but break points
%are smoothened. 
%-> THIS METHOD REQUIRES TO BUILD AN ADJACENCY MATRIX FIRST THAT DEFINES THE
% SPATIAL RELATION BETWEEN NODES, FIRST. 

% Method 2 ('SSCorMeth 2'), uses a fixed increase in the upstream
% elevation, to locally create a very low slope. 

% Method 3: ('SSCorMeth 3'), uses a correlation with the drainage area in order ot derive local slopes. 

%%% INPUTS
% raw_data_sort: Matrix with hydromorphological observation for each reach 
% SSCorMeth: Which correction method to apply
% D: Adjacency matrix (required if SSCOrMeth=1)

%%% Outputs
% raw_data_sort: Matrix with hydromorphological observation for each reach,
% but with CORRECTED SLOPE
% SS_cor_results: Matrix that informs for which reaches slope was
% corrected, and how slope was before and after correction. 


%% 

global ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY % Clear temporary variables

% profile on

switch SSCorMeth
           
            
        case  'SSCorMeth 3'
                       
           ss_val=find(raw_data_sort(:,ID_Slp)>1e-5);% find valid slope indices
           x=raw_data_sort(ss_val,ID_Ad); 
           y=raw_data_sort(ss_val,ID_Slp);

           % Set up fittype and options.
           ft = fittype( 'power1' ); opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); opts.Display = 'Off';
           opts.StartPoint = [0.0528236677872659 -0.169343731781852];

           % Fit model to data.
           [fitresult_slope, gof] = fit( x, y, ft, opts );

       end

%% begin the correction 
SS_cor_results=zeros(length(raw_data_sort),3);
ss=1;
temp=size(raw_data_sort,1)
for sss=1:temp;  %% loop through all reaches 



if raw_data_sort(sss,ID_Slp)<=1e-5 % find reaches with a invalid slope
ss=ss+1;


    switch SSCorMeth
  
    case  'SSCorMeth 1'
         hwb=waitbar(ss/sum(raw_data_sort(:,ID_Slp)<1e-5)); set(hwb,'Name','Correcting Slope'); 
 
         dsEl_sss=raw_data_sort(sss,ID_ElDs); % upstream elevation of current node
         kk=1; % counter for dowsntream reaches  
         dse=nan(1,10000); %storage for downstream reaches
         dse(kk)=sss;
         delta_h=0;
         while delta_h<=0  
            kk=kk+1;
          

            if isempty(find(D(dse(kk-1),:)>0));
                flag=1; break; else  ...  % check if there is a dowsntream node, if not: break
             dse(kk)=find(D(dse(kk-1),:)>0); flag=0; end % find next downstream edge (dse) 
            usel_dse=raw_data_sort(dse(kk),ID_ElUs); % upstream elevation of next downstream edge
            delta_h=(dsEl_sss-usel_dse); %calculate elevation difference
         end 
         
         if flag==0
            newSlope=delta_h./(sum(raw_data_sort(dse(1:kk-1),ID_Length))); %calculate new slope for all required downstream edges 
            raw_data_sort(dse(1:kk-1),ID_Slp)=newSlope  ; %calculate new slope for all required downstream edges
         else 
            raw_data_sort(sss,ID_Slp)=1E-5; %Use very small, random elevation, if there is no downstream node
         end
        
         SS_cor_results(ss,1)=nan;

    case  'SSCorMeth 2'

           %%%% method 2 
           % this is not based in the data, 
           raw_data_sort(sss,ID_Slp)=0.25/raw_data_sort(sss,ID_Length);     
           SS_cor_results(ss,3)=[raw_data_sort(sss,ID_Slp)];
 
    case  'SSCorMeth 3'

           %%%% method 2 
           % this is not based on a extrapolatiioon from Q1.5
           raw_data_sort(sss,ID_Slp)=fitresult_slope(raw_data_sort(sss,ID_Q15));    
           SS_cor_results(ss,3)=[raw_data_sort(sss,ID_Slp)];


    end 
  
end
end

SS_cor_results(SS_cor_results(:,1)==0,:)=[];

close(hwb)
% profile viewer
end

