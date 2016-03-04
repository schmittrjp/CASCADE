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

% tributaries=[632 1348 1398 1425 1769];

Path1_to_out=find(isnan(II(1,:))==0); % find all reaches on the way to the outlet 
% link_dist_from_outlet=fliplr(1:length(Path1_to_out)); % how many links a node is away from the outlet 
link_dist_from_outlet=(1:length(Path1_to_out)); 
ds_path_nodes_full=find(isnan(II(1,:))==0);% Network.Downstream.Path{1,1}{1,3431};
[~,ds_path_nodes_full_sort]=sort(ds_path_nodes_full);
link_dist_from_outlet_sort=link_dist_from_outlet(ds_path_nodes_full_sort); %distance from outlet in number of nodes 
km_dist_from_outlet=Network.Downstream.Distance{1,1}(1,ds_path_nodes_full); %distance from outlet in km 
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

subplot(2,3,1)

kk=0;
  for ll=1:length(AggData)  % plot entire river network 
          hh=line([AggData(ll,ID_FX) AggData(ll,ID_TX)],[AggData(ll,ID_FY) AggData(ll,ID_TY)],'color',[0.7 0.7 0.7],'linewidth',0.1);

  end 
 
 for ll=Path1_to_out % plot the path that is considered 
     kk=kk+1; 
     hold on
   
      hh1=line([AggData(ll,ID_FX) AggData(ll,ID_TX)],[AggData(ll,ID_FY) AggData(ll,ID_TY)],'color','k','linewidth',2); % using a color

 end 
 axis equal square off


                    %    set(c,'yticklabel',num2str(c_class'))
                    %     [fcb,xcb]=ecdf(colors)
                    %     fcb=ceil(fcb*1000)/1000
                    %     cbtick_label=num2str(xcb(find(ismember(fcb,cbtick'))))
                    %     
%                       set(c,'ytick',cbtick,'yticklabel',cbtick_label)
%                     caxis([0 max(colors)])   ;
%                     cbarf(colors,c_class);

% subplot(2,3,2:3)
%     colormap(cmap_2)
%     c=colorbar('horizontal','location','northoutside');
%     set(c,'ytick',[min(c_class2):1/(length(c_class2)):1],'Yticklabel',[num2str(c_class2')],'fontsize',12)
%     axis off
    
temp=Output;
temp(temp==0)=nan;

PP4=(bsxfun(@rdivide, temp', max(temp')))'; % create a new matrix for sediment contribution 
PP4=1-PP4; % how much (cummulative percentage) of the original input in a path is deposited downstream

subplot(2,3,4:6)
%
kk=0;
Path1_to_out=Network.Downstream.Path{1,1}{1,outlet_node_new(1)};
hold on 
h_path=zeros(length(Path1_to_out(1:end-1)),1)

for ptl=Path1_to_out(1:end-1)
        kk=kk+1; 
        [~,c_classkk(kk)]=min((abs(colors(kk)-c_class2)));
   
            ds_path_nodes=Network.Downstream.Path{ptl}{1,outlet_node_new(1)}; % find all downstream nodes in the proper order along the path      
            voidnotes=find(isfinite(Output(ptl,ds_path_nodes))==0); % find which of the downstream nodes are actually not sedimentologically connected
            ds_path_nodes(voidnotes(2:end))=[];
       
       % differentiate between links that reach the outlet and those that
       % do not.. 
             
       if isempty(find(ds_path_nodes==outlet_node_new(1)))==1 

            ds_path_nodes_pos=find(ismember(ds_path_nodes_full,ds_path_nodes)); % find the position of nodes along the entire path  
            xxx=link_dist_from_outlet(ds_path_nodes_pos);
           
%         plot(xxx,[PP3(ptl, ds_path_nodes(1:end-1)) 1],'color',cmap_2(c_classkk(kk),:));
      h_path(ptl)=plot(xxx,[PP4(ptl, ds_path_nodes(1:end-1)) 1],'color',cmap_1(end,:));

        contrib_to_out(ptl,1:4)=[xxx(1) Output(ptl,outlet_node_new(1)) PP4(ptl, ds_path_nodes(end)) c_classkk(kk)]; % Find how much these nodes contribute to the sediemnt output at the basin.outlet 
 
           if kk==258
               klm=1;
           end
        
       else
           
          [~,c_classkk(kk)]=min(abs(PP4(ptl, ds_path_nodes(end-1))-c_class1));
          
          if isnan(c_classkk(kk)); classkk=c_class1(end); end 
           ds_path_nodes(voidnotes(1:end))=[];         % include the first void note, here the cummulative sediment contribution will be 1. See also next lines
           ds_path_nodes_pos=find(ismember(ds_path_nodes_full,ds_path_nodes)); % find the position of nodes along the entire path          
          
           xxx=link_dist_from_outlet(ds_path_nodes_pos);
           
           plot(xxx,[PP4(ptl, ds_path_nodes(1:end))],'color',cmap_1(c_classkk(kk),:));
           contrib_to_out(ptl,1:4)=[xxx(1) Output(ptl,2097) PP4(ptl, ds_path_nodes(end)) c_classkk(kk)]; % Find how much these nodes contribute to the sediemnt output at the basin.outlet 

       end
       
        if any(isnan([PP4(ptl, ds_path_nodes(1:end-1)) 1]));
        klm;
        end
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
 dist_to_out=Network.Downstream.Distance{Path1_to_out};    % distance of all nodes on the paths
 ds_path_nodes=Network.Downstream.Path{1,1}{outlet_node_new(1)};   % sequence of all nodes on the path
    dist_to_out=dist_to_out(ds_path_nodes);                 % reorder distances
    dist_to_out=floor(dist_to_out./1000);

% format the axis    
set(gca,'Xtick',[9:50:259],'Xticklabel',num2str(flipud(dist_to_out(1:50:251)')),'Xlim',[0 kk])   
set(gca,'XAxisLocation','top')
xlabel('Distance from the outlet [km]')
ylabel('Fraction of initial sediment deposited []' )   
% 

%%% overlay a plot with the cumulative sediment budget, or other,
%%% additional information

  
%     %prepare the axis
%     ax1_pos=get(gca,'Position');
%     ax2=axes('Position',ax1_pos); % get axis of main plot 
% 
%     %%%% add a plot with the cummulative sediment balance 
% %             [~,sort_ind]=sort(contrib_to_out(:,1)); % sort the stored values according to the start node 
% %             contrib_to_out_cum=nancumsum(contrib_to_out(sort_ind,2));
% %     h=scatter(ax2,contrib_to_out(sort_ind,1),contrib_to_out_cum,10,'r','filled')
% %     set(ax2,'color','none','Xtick',[],'Xticklabel',[],'Yaxislocation','right','Xlim',[0 kk],'Ycolor','r')
% %     ylabel(ax2,'Cummulative sediment contribution')
% 
%   %%%% add a instead a plot with the d50 along the pathway
% Dmat(isfinite(Tmat_cum)==0)=nan;
% % PP4=bsxfun(@rdivide,Output,nansum(Output));
% PP5=Dmat.*PP;
%     % define which indicator to use for the plot:
%     D50_mat=nanmean(Dmat2); % calculate the D590 using a normal mean 
%     % D50_mat=nanmedian(Dmat); % calculate the D50 using the median
% %     D50_mat=nansum(PP4.*Dmat)./nansum(PP4); % calculate the D50 using a
%     %weighted mean
% k=dii(Path1_to_out);
%     % classifiy k according to grain size 
%       for kkk=1:length(k); [~,k_class(kkk)]=min(abs(k(kkk)-c_class2)); end 
%             
%     % plot 
%        h=scatter(ax2,1:260,zeros(length(k),1),((k_class+1)').^2.5,cmap_2(k_class,:),'filled','linewidth',1,'MarkerEdgeColor','k'); % make the scatter 
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
% %   %%%% add any other variable
% % Dmat(isfinite(Tmat_cum)==0)=nan;
% % PP4=bsxfun(@rdivide,Output,nansum(Output));
% % PP5=Dmat.*PP;
% % D50_mat=nanmean(Dmat); % calculate the D590 using a normal mean 
% % % D50_mat=nanmedian(Dmat); % calculate the D590 using the median
% % %  D50_mat=nansum(PP4.*Dmat)./nansum(PP4); % calculate the D50 using a
% % %weighted mean
% % k=(Ybal_tot(Path1_to_out)); 
% % %  k=(Input(Path1_to_out,2097)); 
% % %  k=PP(Path1_to_out,2095)
% % % k=nancumsum(k);
% % % k = sign(k).*log10(abs(k)); % if required use a sign corrected log-values
% % 
% %     ax1_pos=get(gca,'Position');
% %     ax2=axes('Position',ax1_pos); % get axis of main plot 
%         ylim([-5*10^8 5*10^8])
% % 
% %     ylabel(ax2,'D50 [mm]')
% 
% set(gcf,'Position',get(0,'ScreenSize')*1)
% % 
% % if run_from_top==1
% % cd([folder_id])
% % savefig(h_ds_dep, fig_name)
% % print(gcf,'-depsc2',[fig_name '.eps'],'-r1200')  
% % print(gcf,'-dpdf',[fig_name '.pdf'],'-r1200')    
% % cd(basefolder')
% % end
% figure
% x=dii';
% y=PP4(:,2095)./Distmat_cum(:,2095);
% scatter(x,y)
% ylabel('grain size')
% set(gca,'xscale','log')

%%% use the cursor to find the start node
dcm_obj = datacursormode(gcf);   
datacursormode on

%enable data cursor mode
cursor_info = getCursorInfo(dcm_obj);  %get the properties
set(cursor_info.Target,'LineWidth',2)  %the target handle is Target

c_info = getCursorInfo(dcm_obj);
LineStartNode=Path1_to_out(c_info.Target.XData(1))

% hLine = get(event_obj,'Target')