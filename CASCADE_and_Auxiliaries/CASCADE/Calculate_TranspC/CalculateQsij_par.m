%%   Calculate QSij in each percentile 
% Calculation of intermitency factor for EACH grainsize in EACH reach
% using a parallelized approach

Intermittency_cell=cell(1,length(Dmat)); 
Qs_mean_cell=cell(1,length(Dmat)); 

QSij_cell=cell(1,length(Dmat)); 
   tic 

%%% loop through all pathways 
parfor ii=1:length(AggData) % loop through all observations 

    dii=hydraulicData(ii,2);
%%% loop trough all dowsntream reaches of a pathway 
    rds_ii=[find(isfinite(II(ii,:))==1)]; % reaches downstream of  ii (real indices)
    Qs_mean_cell{ii}=zeros(1,length(Dmat));
    QSij_cell{ii}=zeros(1,length(Dmat));
    Intermittency_cell{ii}=zeros(1,length(Dmat));

        for rds=1:length(rds_ii) % loop through all reaches downstream of reach ii
            jj=rds_ii(rds); %find the real index of reach jj.  
            
            reach_jj=find(AggData(:,ID_FromN)==jj); % position of reach j in the data
           
            if isempty(reach_jj)==0 % check if there are any downstream reaches 

            if numel(reach_jj)>1
            reach_jj=reach_jj(1); %trough the resorting, some reaches are duplicate, use only the first value
            end
                        
%%% loop through the flow percentil
            tcap=zeros(1,length(stats{3,ii})); % tcap stores the trnaport capacit for dii in jj for each flow percentile
            for prc=1:length(stats{3,ii}) 
            
            % hydraulic parameters of the reach, 
                % hydraulic data depend on the current value of the flow percentile
                hjj=stats{5,jj}(prc);
                vjj=stats{6,jj}(prc);
                % morphologic parameters depend on the reach    
                Sjj=AggData(reach_jj,ID_Slp);
                Wjj=AggData(reach_jj,ID_Wac);
            
                if Sjj<10^-4
                    Sjj=10^-4;
                end
                if vjj<10^-5
                    vjj=10^-5;
                end

                t_prc=stats{4,ii}(prc); % number of days within the current percentile
               
                    if Dmat(ii,jj)<=2*10^(-3) % for sand
                        tcap(prc)=Engelund_Hansen(hjj,Dmat(ii,jj),vjj,Wjj,Sjj,9.81)*t_prc; % sand transported  during the flow conditions of t_prc
                    else
                        taucrit=hydraulicData(jj,6);
                        tcap(prc)=Wong_Parker(hjj,Dmat(ii,jj),vjj,Wjj,Sjj,9.81,taucrit)*t_prc; % gravel transported  during the flow conditions of t_prc
                    end   
                                        
                    if isfinite(tcap(prc))==0 || isreal(tcap(prc))==0; tcap(prc)=0; end
           
            end % loop through prc values 
            
            Qs_mean_cell{ii}(jj)=nansum(tcap)./stats{7,jj}; % calculate mean transport capacity [kg/d]

            QSij_cell{ii}(jj)=Qs_mean_cell{ii}(jj).*365; % calculate annual transport capacity [kg/yr]

%%% calculate reference Qs at bankfull stage

 
            t_prc=1; % one day time period 
            hjj=hydraulicData(jj,1); vjj=hydraulicData(jj,4); %use hydraulic parameters at bankfull stage
            if Dmat(ii,jj)<=2*10^(-3) % for sand
                        tcap_1_5=Engelund_Hansen(hjj,Dmat(ii,jj),vjj,Wjj,Sjj,9.81)*t_prc; % Velocity of sediment from i in reach j
            else
                        taucrit=hydraulicData(jj,6);
                        tcap_1_5=Wong_Parker(hjj,Dmat(ii,jj),vjj,Wjj,Sjj,9.81,taucrit)*t_prc; % Velocity of sediment from i in reach j
            end
            
%%% calculate intermittency factor by comparing the annual mean sediment transport with the sedimant trnasport at bankfull stage           
            Intermittency_cell{ii}(jj)=Qs_mean_cell{ii}(jj)./tcap_1_5;
            
            end % if ismepty reach jj
            
 
        end  % loop through all reaches downstream of reach ii 
    
% close(hwb3) % close waitbar
end % loop through all observations 
toc
QSij=(cell2mat(QSij_cell));
QSij=reshape(QSij,size(Dmat))';
QSij(QSij==0)=nan;
meanQSij=nanmean(QSij)./10^3; 