%% plot river network. Connectivity attributes are represented as local color or line-width scale
% % 
% %  if exist('run_from_top')
%     if run_from_top==1
%     h_riv_net=figure('Name',fig_name,'color','w'); 
%     end
% %  else 
% %     run_from_top=0
% %  end
    h_riv_net=figure('Name','River_network','color','w'); 


clear h
clear empty_cat lwdth 

% select which attribute to plot


z=0 ; % counter for categories with entries
%Ybal_tot;
% colordata= t_reachtooutlet; categories=[0:50:1000];
% colordata= Ybal_tot; categories= [-10^9 -10^8 -10^7 -10^6:10^5:10^6 10^9 10^8 10^7]
% colordata=Big_diff_dH; categories= [-1:0.1:1];
% colordata=Big_diff_y; categories=[-8*10^7:10^6:10*10^7];
%         cmap=autumn(sum(categories<=0));
%         cmap=[cmap; gray(1); (parula(sum(categories<0)))];
% colordata=nedge_fastestlink_us; categories=[1:1:20];
%         cmap=autumn(sum(categories<=0));
%         cmap=[cmap; gray(1); (hsv(sum(categories>0)))];
%   colordata= Ybal_tot; categories= [-10^9 -10^8 -10^7 -10^6:10^5:10^6 10^9 10^8 10^7]; cmap=[autumn(sum(categories<=0)); (winter(sum(categories>0)))]; 
% colordata= Ybal_tot3; categories=[ -10^7 10^7 10^12]; cmap=(cbrewer('div','RdYlBu',length(categories)));  

% colordata= That, colordata(isinf(colordata))=nan; categories=[1 2 3 4 5 10 20]; cmap=(cbrewer('div','RdYlBu',length(categories)));  

% colordata= Hyd_out_store(:,7); categories= [10^-3 2*10^-3 0.02 0.065 0.1 0.225]; cmap=(cbrewer('div','BrBG',length(categories))); cmap(1:3,:)=flipud(cmap(1:3,:)) ; 
% colordata= Hyd_out_store(:,1); categories= [0.5 1 2.5 5 10 20]; cmap=(cbrewer('seq','PuBu',length(categories))) ; 
% colordata= Hyd_out_store(:,6); categories= [0.01:0.005:0.1]; cmap=(cbrewer('seq','Reds',length(categories))) ; 

% colordata=ncon_reaches; categories= [10 20 50 100 250 500 750 1000 2000];cmap=[gray(1); (jet(sum(categories>0)))]; lwdth_lim=[0.5 15]; lwdth=lwdth_lim(1):(lwdth_lim(2)-lwdth_lim(1))/length(categories):lwdth_lim(2); output_figname='N_paths_per_reach'
% colordata= t_reachtooutlet; categories=[-9999 5 15 20 30 50 100 200 300 400 500 1000]; cmap=[0 0 0; flipud(cubehelix(sum(categories>0),0.5,-1.5,1,1,[0.29,0.92]))]; output_figname='Connection_time_to_outlet';

% colordata= t_reachtooutlet; categories=[-9999 5 25 50 100 200 300 400 500 1000]; cmap=[0 0 0; (othercolor('Dark23',sum(categories>0)))];...
%     output_figname='Connection_time_to_outlet';
% colordata= t_reachtooutlet; categories=[-9999 5 25 50 100 250 500 1000 2000]; cmap=[0 0 0; (cubehelix(sum(categories>0),0.5,-1.5,1,1,[0.29,0.92]))],lwdth=[0.5 4*ones(1,100)]; output_figname='Connection_time_to_outlet';
%   colordata= Outputs_agg; categories= [10^6  10^7 10^8 10^9 5*10^9 10^10 5*10^10 10^11 5*10^11 10^12 10^13]; cmap=jet(length(categories));  clear lwdth; 
%  colordata=cell2mat(sign_change{3,:}); categories= [0 1 2 3]; cmap=jet(length(categories));  clear lwdth; 

%   colordata= diag(TranspC_mat1); categories= [10^0 10^2 10^3 10^4 10^5 10^6  10^7 10^8 10^9 10^10 10^12]; cmap=(othercolor('BuDOr_18',length(categories))), clear lwdth; 

%   colordata= Tcap2_agg; categories= [10^6  10^7 10^8 10^9 10^10 10^11 10^12]; cmap=jet(length(categories)), clear lwdth; 
% colordata= t_reachtooutlet; categories=[-9999 5 15 20 30 50 100 200 300 400 500 1000 5000]; cmap=flipud(cbrewer('div','RdYlGn',length(categories)+1)); cmap=[0.7 0.7 0.7; cmap]; output_figname='Connection_time_to_outlet';

 colordata=PP4(:,2095); colordata(isnan(colordata))=-9999; categories=[-9999 0.1:0.10:1.00]; cmap=flipud(cbrewer('seq','YlOrBr',length(c_class1)+3)); cmap=[.7 .7 .7; cmap(1:10,:)]; 

if exist('lwdth')==1; yes_linewidth=1; else yes_linewidth=0; end % check if a linewidth attribute is defined 

% if run_from_top==1; clear colordata categories cmap lwdth; load temp_colordata; end % if the script is run from another script delet the colordata and load the colordata form the main script

% [cvals_sort cval_rank]=sort(colordata,'descend');

% the first loop is to find categories that do not contain any data, and to
% adapt the colormap accordingly. 
for nn=1:length(categories)
    
    if nn==1 % find data smaller than the first category, nn=1
        dtp=find(colordata<=categories(nn)); %dtp: data to plot: those data smaller than the current categorie
    else % find data for all other categories nn=2..end
        dtp=find(colordata>categories(nn-1) & colordata<=categories(nn)); %dtp: data to plot: those data larger than the previous, and smaller that the current category 
    end
    
    ptp=find(ismember(raw_data(:,ID_FromN),dtp)); %ptp: position of data (dtp) amongst the original data. this is to get there x and y positions. 
    
    if isempty(ptp)==1
        empty_cat(nn)=1;
    else
        empty_cat(nn)=0;
    end
end 
% if exist('empty_cat')% map data to the entire, initial color-range, even if there are empty categories
% 
%     pos_empty_cat=find(categories) %find which categories are empty 
%     valid_cat=find(empty_cat==0) % find which categories are not empty 
%     
%     cmap(valid_cat(1),:)=cmap(1) % assign the lowest, not-empty-category to the first color value 
%     cmap(valid_cat(end),:)=cmap(end,:) % assign the highest, not-empty-category to the last color value 
%     
%     c_spread=floor(length(categories)/(length(valid_cat)-2))
%     c_spread=valid_cat(1)+1:c_spread:length(categories)-1
%     
%     cmap(valid_cat(2:end-1),:)=cmap(c_spread,:)
%     
%     cmap(valid_cat(1),:)=cmap(1) % assign the lowest, not-empty-category to the first color value 
%     cmap(valid_cat(end),:)=cmap(end,:)  % assign the highest, not-empty-category to the last color value
%     
% end
    % plotting   
for nn=1:length(categories)  
    if nn==1 % find data smaller than the first category, nn=1
        dtp=find(colordata<=categories(nn)); %dtp: data to plot: those data smaller than the current categorie
        else % find data for all other categories nn=2..end
        dtp=find(colordata>categories(nn-1) & colordata<=categories(nn)); %dtp: data to plot: those data smaller than the current categorie
    end
    
    ptp=find(ismember(raw_data(:,ID_FromN),dtp)); %ptp: position of data (dtp) amongst the original data. this is to get there x and y positions. 

      if yes_linewidth==1 % check if an attribute for linewidth is defined 
          
            z=z+1
            for kk=1:length(ptp) % loop through and plot all observations from the current class
              if kk==1  % get the handle form the first line (from the legend)            
              h(z,1)=line([raw_data(ptp(kk),ID_FX) raw_data(ptp(kk),ID_TX)],[raw_data(ptp(kk),ID_FY) raw_data(ptp(kk),ID_TY)],'color',cmap(nn,:),'linewidth',lwdth(nn));
              else 
              line([raw_data(ptp(kk),ID_FX) raw_data(ptp(kk),ID_TX)],[raw_data(ptp(kk),ID_FY) raw_data(ptp(kk),ID_TY)],'color',cmap(nn,:),'linewidth',lwdth(nn));
              end
            end

      else lwdth=raw_data(ptp,10)/1.5;  % if there is no linewidth attribute, scale the lines by Strahler order
       z=z+1
            for kk=1:length(ptp) % loop through and plot all observations from the current class
              if kk==1  % get the handle form the first line (for the legend)            
              h(z,1)=line([raw_data(ptp(kk),ID_FX) raw_data(ptp(kk),ID_TX)],[raw_data(ptp(kk),ID_FY) raw_data(ptp(kk),ID_TY)],'color',cmap(nn,:),'linewidth',lwdth(kk));
              else 
              line([raw_data(ptp(kk),ID_FX) raw_data(ptp(kk),ID_TX)],[raw_data(ptp(kk),ID_FY) raw_data(ptp(kk),ID_TY)],'color',cmap(nn,:),'linewidth',lwdth(kk));
              end
            end
      
      end
      
end
        clear dtp ptp
        axis square equal off 


% prepare the legend 
% create handles for a fake line in order to make a proper legend 
ll_leg=0;
clear h_line
            for ll=find(empty_cat==0)  
                ll_leg=ll_leg+1;
            h_line(ll_leg)=line([raw_data(1,ID_FX) raw_data(1,ID_FX) ],[raw_data(1,ID_FY) raw_data(1,ID_FY)],'color',cmap(ll,:),'linewidth',2);
            end

%         h_orig=get(h,'LineWidth'); % save the original attributes of all line handles 
%         set(h(:),'linewidth',3) % set linewidth to 2pt (so that all lines in the legend are the same)   
        if exist('empty_cat'); cat_for_leg=categories(find(empty_cat==0)  )'; else cat_for_leg=categories'; end  % categories that will go into the legend. 
        
        if max(cat_for_leg)>10^7; leg_ent=num2str(cat_for_leg,'%10.2E'); else leg_ent=num2str(cat_for_leg); end; % prepare the text for the legend
        hleg=legend(h_line,leg_ent); % create the legend 
        set(hleg,'FontSize',14); % set larger fontsize 

set(0,'DefaultFigureColormap',cmap);
cb=colorbar; 
caxis([min(colordata') max(colordata')])
cstep=(max(colordata')-min(colordata'))/length(categories);
% set(cb,'Ytick',min(colordata'):cstep:max(categories),'Yticklabel',num2str([min(colordata'); cat_for_leg; max(colordata')]));
set(cb,'Ytick',[min(colordata'):cstep:max(colordata') ],'Yticklabel',num2str([min(colordata'); cat_for_leg; max(colordata')],'%2E'));

set(gcf,'Position',get(0,'ScreenSize'))

if run_from_top==1
cd([folder_id])
savefig(h_riv_net, fig_name)
print(gcf,'-depsc2',[fig_name '.eps'],'-r1200')  
print(gcf,'-dpdf',[fig_name '.pdf'],'-r1200')    
cd(basefolder')
end