%% plot the bottlenecks according to their type: Due to competition, or local insufficient transport capacity
% figure 
%
n_path_min=5; % minimum number of path interrupted in a reach

% PLOT_River_network % plot river network with an additional attribute

%       for ll=length(AggData)  % plot river network in black
%               hh=line([AggData(ll,ID_FX) AggData(ll,ID_TX)],[AggData(ll,ID_FY) AggData(ll,ID_TY)],'color',[0.3 0.3 0.3],'linewidth',AggData(ll,10)/1.5);
% 
%       end 
%       axis square equal off 
%

Discon_scatter{1}=zeros(size(Discon_for_mbal,1),4);
Discon_scatter{2}=zeros(size(Discon_for_mbal,1),4);

uni_mbal=unique(Discon_for_mbal(:,2)); % unique 

if isempty(uni_mbal)==0 % check for interruptions due to competition

        for klm=1:length(uni_mbal)
           uni_mbal(klm,2)=sum((Discon_for_mbal(:,2)==uni_mbal(klm))); 
        end
        
[~, temp]=sort(uni_mbal(:,2),'descend'); 
uni_mbal=uni_mbal(temp,:); 
% Discon_scatter{2}=uni_mbal(1:n_path_for_plot,1:2);
Discon_scatter{1}=uni_mbal(1:end,1:2);

Discon_scatter{1}(Discon_scatter{1}(:,2)<n_path_min,:)=[]

clear temp        

kk(1)=1 % for plotting       
else
kk(1)=[]    

end         
 
uni_local=unique(Discon_for_local(:,2)); 
   uni_local(uni_local==0)=[] % check and remove zero values from the discon scatter 


if isempty(uni_local)==0 % check for diconnections due to local energy limit
        for klm=1:length(uni_local)
           uni_local(klm,2)=sum((Discon_for_local(:,2)==uni_local(klm))); 

        end
        
[~, temp]=sort(uni_local(:,2),'descend'); 
uni_local=uni_local(temp,:); 
Discon_scatter{2}(1:length(uni_local),1:2)=uni_local(1:end,1:2); 

Discon_scatter{2}(Discon_scatter{2}(:,2)<n_path_min,:)=[]
kk(2)=2
else


end

% clear Discon_scatter

 for kkk=kk % loop trhough local and massbalnce interrupted paths 
    for discon=1:size(Discon_scatter{kkk},1) % loop through disconnected due to the mass balance 
     node=find(AggData(:,ID_FromN)==Discon_scatter{1,kkk}(discon,1));
     Discon_scatter{kkk}(discon,3)= AggData(node,ID_FX);    % find X coord 
     Discon_scatter{kkk}(discon,4)= AggData(node,ID_FY) ;   % find Y coord
    end
 end
hold on 

if kk(1)==1 % Round mrkers indicate disconectivity for local competition 
 scatter(Discon_scatter{1}(:,3),Discon_scatter{1}(:,4),Discon_scatter{1}(:,2)*5,'b','filled')
end
hold on

if length(kk)==2 % square markers indicate disconnectivity for local energy limit
 scatter(Discon_scatter{2}(:,3),Discon_scatter{2}(:,4),Discon_scatter{2}(:,2)*5,'r','s','filled')
end

%% make some fake points for the legend
scatterPointSize=[5 15 30]; % which dot sizes to create 

scatterx=(0.8:0.1:1)*max(Discon_scatter{kkk}(discon,3)); % x position of dots
scatter1y=ones(size(scatterPointSize))*max(raw_data(:,ID_FY)); % y position for dots (interruption for mass balance)
scatter2y=ones(size(scatterPointSize))*max(raw_data(:,ID_FY))*0.95; % y position for dots (interruptions for local effects). 

hold on 

scatter(scatterx,scatter1y,scatterPointSize*5,'b','filled')
scatter(scatterx,scatter2y,scatterPointSize*5,'r','s','filled')







 