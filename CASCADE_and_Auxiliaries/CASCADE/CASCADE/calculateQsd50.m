%%    3: Calculate QSd50 in each reach 
% Calculation of intermitency factor for EACH grainsize in EACH reach

Intermittency_mat=zeros(size(Dmat)); 
Qsd50=nan(1,size(Dmat,2)); 
Taucrit_mat=nan(size(Dmat));

%%% loop through all reaches
for iii=1:length(d50) % loop through all observations 
  iii;
  if iii==2139; 
      asdsad=1; 
  end 
%  multiWaitbar('Calculating percentile sediment transport' ,ii/length(AggData));

    dii=d50(iii);
%%% loop trough all dowsntream reaches of a pathway 
           
            if isempty(reach_jj)==0 % check if there are any downstream reaches 

            if numel(reach_jj)>1
            reach_jj=reach_jj(1); %trough the resorting, some reaches are duplicate, use only the first value
            end
            
            clear tcap 
            
%%% loop through the flow percentile 
            for prc=1:length(stats{3,iii}) 
             qsd50=zeros(1,length(stats{3,iii}));
            % hydraulic parameters of the reach, 
                % hydraulic data depend on the current value of the flow percentile
                hii=stats{5,iii}(prc);
                vii=stats{6,iii}(prc);
                % morphologic parameters dependend on the downstream reach    
                Sii=AggData(iii,ID_Slp);
                Wii=AggData(iii,ID_Wac);
                               
                % calculate critical shear stress for dii in downstram reach j     
                    Rh=(hii*Wii)/(2*hii+Wii);

                    if Rh/dii<10  % for low water levels in comparison to the grain size, use zthe formula proposed by Suszka L (1991) Modification of transport rate formula for steep channels. Fluvial Hydraulics of Mountain Regions, Lecture Notes in Earth Sciences., eds Armanini PA, Silvio PGD (Springer Berlin Heidelberg), pp 59–70. Available at: http://link.springer.com/chapter/10.1007/BFb0011182 [Accessed March 2, 2015].
                        taucrit=0.0851*(hii/dii)^(-0.0261);
                    else  
                        Re_dii=(1.6*9.81*dii)^0.5*dii/10^(-6); % particle Reynolds number
                        taucrit=(0.22*Re_dii^(-0.6)+0.06*10^(-7.7*Re_dii^(-0.6)));
                    end
%                     Taucrit_mat(ii,jj)=taucrit; 
                
                if Sii<10^-4
                    Sii=10^-4;
                end
                if vii<10^-5
                    vii=10^-5;
                end
                
%%% calculate sediment transport capacity for grain size di in reach j for flow percentile n                  
                    if dii<=2*10^(-3) % for sand
                        qsd50(prc)=Engelund_Hansen(hii,dii,vii,Wii,Sii,9.81); % sand transported  during the flow conditions of t_prc
                    else
%                          taucrit=hydraulicData(jj,6);
                        qsd50(prc)=Wong_Parker(hii,dii,vii,Wii,Sii,9.81,taucrit); % gravel transported  during the flow conditions of t_prc
                    end   
                                        
                    if isfinite(qsd50(prc))==0 || isreal(qsd50(prc))==0; qsd50(prc)=0; end              

             qsd50(prc)=qsd50(prc)*stats{4,iii}(prc); % kg total transport over all  all observation in the current percentile (in kg)  
           
            end % loop through prc values 
      
            Qsd50(iii)=sum(qsd50)./sum(stats{4,iii}); %divide by total number of days on record [kg/d]
            Qsd50(iii)=Qsd50(iii)*365; % scale to kg/yr
            
        end  % loop through all reaches downstream of reach ii 
    
end % loop through all observations 

% multiWaitbar('Closeall');
    


    
