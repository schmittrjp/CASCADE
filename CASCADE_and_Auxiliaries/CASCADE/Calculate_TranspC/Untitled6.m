%%    3: Calculate QSij in each percentile 
% Calculation of intermitency factor for EACH grainsize in EACH reach

profile on 
Intermittency_mat=zeros(size(Dmat)); 
Qsij=nan(size(Dmat)); 
Taucrit_mat=nan(size(Dmat));

%%% loop through all pathways 
for ii=1:length(AggData) % loop through all observations 
  ii; 
 multiWaitbar('Calculating percentile sediment transport' ,ii/length(AggData));

    dii=hydraulicData(ii,2);
%%% loop trough all dowsntream reaches of a pathway 
        for jj=find(isfinite(II(ii,:))==1) % loop through all reaches downstream of reach ii
            reach_jj=find(AggData(:,ID_FromN)==jj); % position of reach j in the data
           
            if isempty(reach_jj)==0 % check if there are any downstream reaches 

            if numel(reach_jj)>1
            reach_jj=reach_jj(1); %trough the resorting, some reaches are duplicate, use only the first value
            end
            
            clear tcap
            qsij=zeros(1,length(stats{3,ii}));
%%% loop through the flow percentile 
            for prc=1:length(stats{3,ii}) 
                
                
                if prc==5
                    ads=21; 
                end 
            
            % hydraulic parameters of the reach, 
                % hydraulic data depend on the current value of the flow percentile
                hjj=stats{5,jj}(prc);
                vjj=stats{6,jj}(prc);
                % morphologic parameters dependend on the downstream reach    
                Sjj=AggData(reach_jj,ID_Slp);
                Wjj=AggData(reach_jj,ID_Wac);
                               
                % calculate critical shear stress for dii in downstram reach j     
                    Rh=(hjj*Wjj)/(2*hjj+Wjj);

                    if Rh/dii<10  % for low water levels in comparison to the grain size, use zthe formula proposed by Suszka L (1991) Modification of transport rate formula for steep channels. Fluvial Hydraulics of Mountain Regions, Lecture Notes in Earth Sciences., eds Armanini PA, Silvio PGD (Springer Berlin Heidelberg), pp 59–70. Available at: http://link.springer.com/chapter/10.1007/BFb0011182 [Accessed March 2, 2015].
                        taucrit=0.0851*(hjj/dii)^(-0.0261);
                    else  
                        Re_dii=(1.6*9.81*dii)^0.5*dii/10^(-6); % particle Reynolds number
                        taucrit=(0.22*Re_dii^(-0.6)+0.06*10^(-7.7*Re_dii^(-0.6)));
                    end
                    Taucrit_mat(ii,jj)=taucrit; 

                
                if Sjj<10^-4
                    Sjj=10^-4;
                end
                if vjj<10^-5
                    vjj=10^-5;
                end
                
%%% calculate sediment transport capacity for grain size di in reach j for flow percentile n                  
                    if dii<=2*10^(-3) % for sand
                        qsij(prc)=Engelund_Hansen(hjj,dii,vjj,Wjj,Sjj,9.81); % sand transported  during the flow conditions of t_prc
                    else
%                          taucrit=hydraulicData(jj,6);
                        qsij(prc)=Wong_Parker(hjj,dii,vjj,Wjj,Sjj,9.81,taucrit); % gravel transported  during the flow conditions of t_prc
                    end   
                                        
                    if isfinite(qsij(prc))==0 || isreal(qsij(prc))==0; qsij(prc)=0; end              
%              global
%              [qsij(prc)] = calcQsij(dii,hjj,vjj,Sjj,Wjj,taucrit,1); % in kg/d 

             qsij(prc)=qsij(prc)*stats{4,ii}(prc); % kg total transport over all  all observation in the current percentile (in kg)  
%              
            end % loop through prc values 
      
            Qsij(ii,jj)=sum(qsij)./sum(stats{4,ii}); %divide by total number of days on record [kg/d]
            Qsij(ii,jj)=Qsij(ii,jj)*365; % scale to kg/yr
            
%%% calculate reference Qs at bankfull stage

            h1_5=hydraulicData(jj,1); v1_5=hydraulicData(jj,4); %use hydraulic parameters at bankfull stage
            Qsij1_5=calcQsij(dii,h1_5,v1_5,Sjj,Wjj,taucrit,365);
 
%%% calculate intermittency factor by comparing the annual mean sediment transport with the sedimant trnasport at bankfull stage           
            Intermittency_mat(ii,jj)=Qsij(ii,jj)./Qsij1_5;
            
            if Qsij(ii,jj)==0 ;          
                break  % if the transport capacity for di in any node j is 0, then there is no need to calculate the transportcapacity in further downstream nodes. 
            end 
            
            end % if isempty reach jj
            
 
        end  % loop through all reaches downstream of reach ii 
    
% close(hwb3) % close waitbar
end % loop through all observations 

multiWaitbar('Closeall');
    
profile viewer

Qsij(Qsij==0)=nan;
  
  %%
meanQSij=nanmean(Qsij)./10^3; 
 
 QS_classes=([ 1E4 1E5 1E6 1E7 1E8 1E9 1E10]);

 networkPlotter(AggData,meanQSij,QS_classes, AggData(:,ID_StrO),'Q_S Classes','Q_S [mio t / yr ]')

    
