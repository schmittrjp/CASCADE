%% General aim is to show how sediment is deposited along the various pathways from node 1 to the outlet. 
% therefor, plot first a map with the path from 1-to-outlet and than the
% plot with the deposition of each path, and  the total sediment transport 
%  if exist('run_from_top') 
%     h_ds_dep=figure('Name',fig_name,'color','w'); 
%  else 
%     run_from_top=0;
%  end
figure('Name','Sed contribution to outlet scenario','color','w')
clear contrib_to_out

% define confluences with major tributaries

tributaries=[632 1120 1398 1425 1769];

Path1_to_out=Network.Downstream.Path{1,1}{1,outlet_node_new(1)}; % find all reaches on the way to the outlet 
% link_dist_from_outlet=fliplr(1:length(Path1_to_out)); % how many links a node is away from the outlet 
temp=Network.Downstream.Path{1,1}(Path1_to_out);
link_dist_from_outlet=cellfun('length',temp); clear temp 
% Path1_to_out=Network.Downstream.Path{1,1}{1,3431};
[~,Path1_to_out_sort]=sort(Path1_to_out);
link_dist_from_outlet_sort=link_dist_from_outlet(Path1_to_out_sort); %distance from outlet in number of nodes 
km_dist_from_outlet=Network.Downstream.Distance{1,1}(1,Path1_to_out); %distance from outlet in km 
% color definition
colors=Dmat(1,Path1_to_out);
% [d_sort,value_sort]=sort(colors,'descend');
c_class1=[0.1:0.1:1];
c_class2=[10^-3 2*10^-3 0.02 0.065 0.1 0.225];

cmap_1=flipud(cbrewer('seq','YlOrBr',length(c_class1)+3));%cmap1 can be used to color the deposition lines by the final value of deposition
%     cmap_1=cmap_1(1:10,:);
    cmap_1=[cmap_1(1:10,:); [0.7 0.7 0.7]];
% cmap_2=(parula(length(c_class2))); % cmap 2 colors deposition lines by grain size 
 cmap_2=flipud((cbrewer('seq','Greys',length(c_class2))));
set(0,'DefaultFigureColormap',cmap_1);

axis equal square off
    
temp=Output;
temp(temp==0)=nan;

PP4=(bsxfun(@rdivide, Output_Scn3', nanmax(Output_Scn3')))';
PP4_1=(bsxfun(@rdivide, Output', nanmax(Output')))';
% create a new matrix for sediment contribution 
PP4=1-PP4; % how much (cummulative percentage) of the original input in a path is deposited downstream

line_ax=subplot(2,3,4:6);
%
kk=0;
Path1_to_out=Network.Downstream.Path{1,1}{1,outlet_node_new(1)};
hold on 

for ptl=(Path1_to_out(1:end-1))
        kk=kk+1; 
                   if kk==26
               klm=1;
           end
        [~,c_classkk(kk)]=min((abs(colors(kk)-c_class2)));
   
            ds_path_nodes=Network.Downstream.Path{ptl,1}{1,outlet_node_new(1)}; % find all downstream nodes in the proper order along the path      
            voidnotes=find(isnan(Output(ptl,ds_path_nodes)) | Output(ptl,ds_path_nodes)==0); % find which of the downstream nodes are actually not sedimentologically connected
            ds_path_nodes(voidnotes(2:end))=[];
       
       % differentiate between links that reach the outlet and those that
       % do not.. 
             
       if isempty(find(ds_path_nodes==outlet_node_new(1)))==1 

            ds_path_nodes_pos=find(ismember(Path1_to_out,ds_path_nodes)); % find the position of nodes along the entire path  
            xxx=link_dist_from_outlet(ds_path_nodes_pos);
           
            if any(diff(xxx)<0)==1
            asdsad=0; 
            end
            
%         plot(xxx,[PP3(ptl, ds_path_nodes(1:end-1)) 1],'color',cmap_2(c_classkk(kk),:));
          h(kk)=plot(xxx,[PP4(ptl, ds_path_nodes(1:end-1)) 1],'color',cmap_1(end,:));

        contrib_to_out(ptl,1:5)=[xxx(1) Output(ptl,outlet_node_new(1)) PP4(ptl, ds_path_nodes(end)) c_classkk(kk) ptl]; % Find how much these nodes contribute to the sediemnt output at the basin outlet 
 

        
       else
           
          [~,c_classkk(kk)]=min(abs(PP4(ptl, ds_path_nodes(end-1))-c_class1));
          
          if isnan(c_classkk(kk)); classkk=c_class1(end); end 
%            ds_path_nodes(voidnotes(1:end))=[];         % include the first void note, here the cummulative sediment contribution will be 1. See also next lines
           ds_path_nodes_pos=find(ismember(Path1_to_out,ds_path_nodes)); % find the position of nodes along the entire path          
          
           xxx=link_dist_from_outlet(ds_path_nodes_pos);
           
           h(kk)=plot(xxx,[PP4(ptl, ds_path_nodes(1:end))],'color',cmap_1(c_classkk(kk),:));
           contrib_to_out(ptl,1:4)=[xxx(1) Output(ptl,outlet_node_new(1)) PP4(ptl, ds_path_nodes(end)) c_classkk(kk)]; % Find how much these nodes contribute to the sediemnt output at the basin.outlet 

       end
       
        if any(isnan([PP4(ptl, ds_path_nodes(1:end-1)) 1]));
        klm;
        end
        
   % plot markers for tributaries     
   if sum(ptl==(tributaries))==1
   patch([xxx(1)-2,xxx(1),xxx(1)+2],[0.2 0 0.2],'k')
   end 
        
end


contrib_to_out((contrib_to_out(:,1)==0 & contrib_to_out(:,2)==0),: )=[];
set(0,'DefaultFigureColormap',cmap_1);

xlim([0 kk])  
% 
 % GET the distances of the nodes on the plotted path.
 % this is to label the x axis with distances, not only with node numbers 
 dist_from_out=Network.Downstream.Distance{Path1_to_out,1};    % distance of all nodes on the paths
 ds_path_nodes=Network.Downstream.Path{1,1}{1,outlet_node_new(1)};   % sequence of all nodes on the path
    dist_from_out=dist_from_out(ds_path_nodes);                 % reorder distances
    dist_from_out=floor(dist_from_out./1000);

% format the axis    
xticks=[0:floor(length(ds_path_nodes)/10):length(ds_path_nodes)];
set(gca,'Xtick',xticks,'Xticklabel',num2str(dist_from_out([1 xticks(2:end)])'),'Xlim',[0 kk])   
set(gca,'XAxisLocation','top')
xlabel('Distance from the outlet [km]')
ylabel('Fraction of initial sediment deposited []' )   
% 

% %% make a second set of axes to display source grain size
%   %prepare the axis
%     ax1_pos=get(gca,'Position');
%     ax2=axes('Position',ax1_pos); % get axis of main plot 
% 
% k=Dmat(Path1_to_out,outlet_node_new(1));
%     % classifiy k according to grain size 
%       for kkk=1:length(k); [~,k_class(kkk)]=min(abs(k(kkk)-c_class2)); end 
%             
%     % plot 
%        h=scatter(ax2,1:length(k),zeros(length(k),1),((k_class+1)').^2.5,cmap_2(k_class,:),'filled','linewidth',1,'MarkerEdgeColor','k'); % make the scatter 
%                     % make legend, loop thourh all used marker classes and
%                     % create a handle for each marker class
%                         hold on; clear h_kk
%                         kkkk_ind=0; 
%                         for kkkk=(unique(k_class))
%                             kkkk_ind=kkkk_ind+1; 
%                         h_kk(kkkk_ind)=scatter(ax2,-1000,k_class(kkkk),(k_class(kkkk)+1).^3,cmap_2((kkkk),:),'filled','linewidth',1,'MarkerEdgeColor','k'); 
%                         end
%                         hold off
%                         kk_legend_text=num2str((c_class2(unique(k_class))*1000)');
%                         h_leg=legend(h_kk,{kk_legend_text},'FontSize',12,'Position',[0.8,0.75,0,0]);
% %                         legendTitle(h_leg,'Grain size [mm]','Fontsize',12,'Fontweight','bold')
%        set(ax2,'color','none','Xtick',[],'Xticklabel',[],'Ytick',[],'Yticklabel',[],'Yaxislocation','right','Xlim',[0 kk],'Ylim',[0 1]);%,'Ycolor','r');%,'Yscale','log')    ylim([-5*10^11 5*10^11])
%    
% %               set(ax2,'color','none','Xtick',[],'Xticklabel',[],'Ytick',[],'Yticklabel',[],'Yaxislocation','right','Xlim',[0 50],'Ylim',[0 1]);%,'Ycolor','r');%,'Yscale','log')    ylim([-5*10^11 5*10^11])
% 
%     % Make colorbar, use a different axis, as not to disturb the other two
%     % plots 
%     cax_pos=ax1_pos+[ax1_pos(3)+0.015, 0, -0.97*ax1_pos(3), 0.1*ax1_pos(4)]; % position of the new colorbar axis
%     ax3=axes('Position',cax_pos);
%     set(ax3,'color','none','Ytick',[],'yticklabel',[],'Xtick',[],'Xticklabel',[],'Fontsize',12)
%     cb=colorbar(ax3,'Position',cax_pos);
%     caxis([0 1.1]);
%     set(cb,'Ytick',[0:0.1:1.0 1.05],'yticklabel',num2str(flipud([0:0.1:1]')))  
%     ylabel(cb,['Fraction of sediment' '\newline' 'delivered to outlet []'],'Fontsize',12) 
%  
%     set(findall(gcf,'type','text'),'FontSize',12)
% 
%     
%     ylim([-5*10^8 5*10^8])
% % 
% %     ylabel(ax2,'D50 [mm]')
% 
% set(gcf,'Position',get(0,'ScreenSize')*1)

%% enable data cursor mode: clicking on a line wil show its handle
%%% use the cursor to find the start node
global h_selected
set(h, 'ButtonDownFcn', {@myLineCallback, h, Input, Output, PP4, Network, outlet_node_new, TranspC_mat2 });
h_selected