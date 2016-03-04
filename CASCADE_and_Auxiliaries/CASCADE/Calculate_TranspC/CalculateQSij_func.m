function [ QSij  Intermittency] = CalculateQSij_func(ds, i, j, Sjj, Wjj, taucrit, stats, vbf, hbf)

% CalculateQSij_func:  Calculate QSij and intermittency factor for a given grain size ds in a single or multiple reaches jj. 
%   dii=hydraulicData(ii,2); 
%   Intermittency_cell=cell(1,length(Dmat)); 
%   Qs_mean_cell=cell(1,length(Dmat)); 
%   QSij_cell=cell(1,length(Dmat)); 

% Inputs: 
% ds: grain size of source s 
% i: reach in which source s is located
% j: reach or reaches downstream of s
% Wjj: Active channel width of reach(es) jj
% Sjj: Slope of reach(es) jj
% taucritjj: critical shear stress in reach(es) jj
% stats: Hydraulic conditions in reach(es) jj
% vbf, hbf: Bankfull velocity and water level 

% Outputs: 

%%% Initialize storage
QSij=zeros(1,length(j));
Qs_mean=zeros(1,length(j));
Intermittency=zeros(1,length(j));
 
%%% loop trough all reaches jj

        for jj=1:length(j) % loop through all reaches j
            %find the real index of reach jj.  
%%% loop through the flow percentiles in reach jj
            tcap=zeros(1,length(stats{3,jj})); % tcap stores the trnaport capacit for dii in jj for each flow percentile
            for prc=1:length(stats{3,jj}) 
            
            % hydraulic parameters of reach jj, 
                hjj=stats{5,jj}(prc);
                vjj=stats{6,jj}(prc);
                            
                if Sjj<10^-4
                    Sjj=10^-4;
                end
                if vjj<10^-5
                    vjj=10^-5;
                end

                t_prc=stats{4,jj}(prc); % number of days within the current percentile
               
                    if ds<=2*10^(-3) % for sand
                        tcap(prc)=Engelund_Hansen(hjj,ds,vjj,Wjj(jj),Sjj(jj),9.81)*t_prc; % sand transported  during the flow conditions of t_prc
                    else
                        tcap(prc)=Wong_Parker(hjj,ds,vjj,Wjj(jj),Sjj(jj),9.81,taucrit(jj))*t_prc; % gravel transported  during the flow conditions of t_prc
                    end   
                                        
                    if isfinite(tcap(prc))==0 || isreal(tcap(prc))==0; tcap(prc)=0; end
           
            end % loop through prc values 
            
            Qs_mean(jj)=nansum(tcap)./stats{7,jj}; % calculate mean transport capacity [kg/d]
            QSij(jj)=Qs_mean(jj).*365; % calculate annual transport capacity [kg/yr]

%%% calculate reference Qs at bankfull stage
 
            t_prc=1; % one day time period 
            %use hydraulic parameters at bankfull stage
            if ds<=2*10^(-3) % for sand
                        tcap_1_5=Engelund_Hansen(hbf(jj),ds,vbf(jj),Wjj(jj),Sjj(jj),9.81)*t_prc; % Velocity of sediment from i in reach j
            else
                        tcap_1_5=Wong_Parker(hbf(jj),ds,vbf(jj),Wjj(jj),Sjj(jj),9.81,taucrit(jj))*t_prc; % Velocity of sediment from i in reach j
            end
            
%%% calculate intermittency factor by comparing the annual mean sediment transport with the sedimant trnasport at bankfull stage           
            Intermittency(jj)=Qs_mean(jj)./tcap_1_5;
  
        end  % loop through all reaches downstream of source s 

end

