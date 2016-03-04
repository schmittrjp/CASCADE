

clear Agg_data sectMS
% % find all  headwater nodes (definition: a headwater nodde is a node that is only a from-node, but no to-node)
% 
% disp('Identifying sections')
% 
% headwaters=raw_data_sort_ex(~ismember(raw_data_sort_ex(:,ID_FromN),raw_data_sort_ex(:,ID_ToN)),ID_FromN); % very rapid approach from: http://www.mathworks.com/matlabcentral/newsreader/view_thread/300928
% 
% % find all confluences (definition: a confluence node is a node that is to-node for multiple edges.) 
% 
% [temp, ~]=histcounts(raw_data_sort_ex(:,ID_ToN),1:length(raw_data_sort_ex));
% confluences=raw_data_sort_ex(temp>1,ID_FromN); clear temp
% 
% sect_start_nodes=[headwaters; confluences];
% 
% % write section membership table: which nodes belong to which section? 
% 
% 
% for iii=1:length(sect_start_nodes)
%     sectMS{iii,1}=sect_start_nodes(iii);
% end 
% 
%     h = waitbar(0,'Initializing waitbar...');
% 
% 
% for ll=1:length(sect_start_nodes)
%    
%     waitbar(ll/length(sect_start_nodes),h,'Deriviing sections')
% 
%    [~, sectMS{ll,2}, ]=graphshortestpath(D,sect_start_nodes(ll),outlet_node); % find  the path from each section start node to the outlet.
%    temp=find(ismember(sectMS{ll,2},sect_start_nodes),2); % along that path find which nodes mark the beginning of another section
%   
%    if length(temp) > 1 % check if there is any downstream confluence (not the case, e.g., for reaches close to the basin outlet)
%    next_confluence=temp(2); % the second value in temp is the next confluence, the first value is actually the current node
%  
%    temp2=sectMS{ll,2}; % get the entire path 
%    temp2(next_confluence:end)=[]; % delete all nodes downstream of thenext confluence 
%    sectMS{ll,2}=temp2; % save the resulting section
%    end  
%    
% end

%% Aggregation, identification of subsections
%loop through all sections and aggregate reaches that belong to the same HYMO type 

% this step requires first to define also connectivity in the upstream
% direction. 
[DUs,~]=write_adj_matrix( raw_data_sort_ex(:,ID_ToN),raw_data_sort_ex(:,ID_FromN),raw_data_sort_ex(:,ID_Length));

%%

disp('Aggregating reaches based on Hymo class')

Agg_data=zeros(length(raw_data_sort_ex),22);
dissolve=zeros(length(raw_data_sort_ex),1);
rand_col=zeros(length(raw_data_sort_ex),1);

h = waitbar(0);
i=0;
idSubSect=0; % counter for the subsections. 
for s=1:length(sectMS)

    waitbar(s/length(sectMS),h,'Identifying sections')

    
     sectionMem=sectMS{s,2}; % all members of a section 

     % The measurements of HYMO parameters can be noisy, if reach length is
     % short. Therefor it is best to apply an moving average filter to filter
     % measurements. Using a filter with a filter window of length
     % n=movavgWindow, means loosing observations for n-1 nodes. Therefor, the
     % section need to be increased, to include n-1 nodes upstream of the
     % section's first, and n-1 sections after the section's last node. 
     
     % define the length of the moving average window    
     movavgWindow=floor(length(sectMS{s,2})/8); % this is a parameter for further experiments. 
     if movavgWindow==0; movavgWindow=1; end; 
     
     % find n-1 upstream nodes    
     sectionUsNodes=findNUsNodes(DUs,sectionMem(1),movavgWindow-1,raw_data_sort_ex(:,ID_Ad));
     
     % find n-1 downstream nodes    
     sectionDsNodes=findNDsNodes(D,sectionMem(end),movavgWindow-1,raw_data_sort_ex(:,ID_Ad));

     % Enlarge the section by these nodes. 
     sectionMem=[sectionUsNodes' sectionMem sectionDsNodes'];
     
     % if there are no upstream nodes, delete the nans resulting from the
     % previous step
     sectionMem(isnan(sectionMem))=[];
     

    % detect downstream change in parameters 
    AbsChangeSlp=raw_data_sort_ex(sectionMem,ID_Slp);
    AbsChangeFldPlnWdth=raw_data_sort_ex(sectionMem,ID_FldPlnWdth);
    AbsChangeConf=raw_data_sort_ex(sectionMem,ID_Conf);
    
   %%
    if length(sectMS{s,2})>1
        
        % prepare input data for k-means clustering for the members of the
        % current section. 
        clust_data=[AbsChangeSlp  AbsChangeFldPlnWdth AbsChangeConf];

        % normalize input data between 0 and 1 
        clust_data_norm=bsxfun(@rdivide, clust_data,max(clust_data));
        
        % calculate the relative downstream change in each parameter. 
        ds_clust_data_change=abs(diff(clust_data_norm));
        
        % sum the normalized change for the parameters -> measure for how
        % strongly HYMO parameters change between reaches. 
        ds_clust_data_change_tot=[0; sum(abs(diff(clust_data_norm)),2)];
        
        clust_data_norm=[clust_data_norm ds_clust_data_change_tot];
        
        % the data are quite noisy because of the small length of
        % the reaches. Use a moving average to smooth the data 
        
        clust_data_norm_mov=tsmovavg(clust_data_norm,'s',movavgWindow,1); 
       
        % delete section members for which there is no moving average value
        % (i.e., nodes at the beginning of a section that has no upstream
        % nodes). 
        sectionMem(isnan(clust_data_norm_mov(:,1)))=[];
        clust_data_norm_mov(isnan(clust_data_norm_mov(:,1)),:)=[];
        
        % delete section members, that where added dowsntream after the end
        % of the original section. Length of this deletion is movavgWindow-2
        sectionMem(end-(movavgWindow-2):end)=[];
        clust_data_norm_mov(end-(movavgWindow-2):end,:)=[];
        
        
     % clustering 
     %         if length(sectMS{s,2})<10
%         n_clust=2;
%         else
%         n_clust=ceil(length(sectMS{s,2})/10);
%         end  

        if length(sectMS{s,2})<5
         n_clust=1;
        else
          [C,~]= subclust(clust_data_norm_mov,1);
         n_clust=size(C,1);
         end
        
        
        kCid=kmeans(clust_data_norm_mov,n_clust);

%         test=AbsChangeSlp(1:end-1)./AbsChangeSlp(2:end);

        type_change=abs(diff(kCid)); 
        temp=type_change>0;
        type_change=[1 temp']; clear temp; 
        
    else 
     type_change=1;   
    end
    
    % every group of nodes belonging to the same HYMO class
    % gets its own subsection
   
    subsectionID=cumsum(type_change);
    
%         figure
%             subplot(2,1,1)
%             stairs(clust_data_norm_mov)
%             legend('slp','fldpln','conf','chnage')
%             subplot(2,1,2)
%             stairs(kCid);ylim([0 n_clust+1])
   %%     
   
    %
    for ss=unique(subsectionID)
        i=i+1;idSubSect=idSubSect+1;
        
    subsectMem=sectionMem(subsectionID==ss);
    
    Agg_data(i,ID_arcid)=raw_data_sort_ex(subsectMem(1),ID_arcid);
    Agg_data(i,ID_FromN)=raw_data_sort_ex(subsectMem(1),ID_FromN);
    Agg_data(i,ID_ToN)=raw_data_sort_ex(subsectMem(end),ID_ToN);
    Agg_data(i,ID_ElUs)=raw_data_sort_ex(subsectMem(1),ID_ElUs);
    Agg_data(i,ID_ElUsRaw)=raw_data_sort_ex(subsectMem(1),ID_ElUsRaw);
    Agg_data(i,ID_ElDs)=raw_data_sort_ex(subsectMem(end),ID_ElDs);
    Agg_data(i,ID_ElDsRaw)=raw_data_sort_ex(subsectMem(end),ID_ElDsRaw);
    Agg_data(i,ID_Slp)=mean(raw_data_sort_ex(subsectMem,ID_Slp));
    Agg_data(i,ID_SlpRaw)=mean(raw_data_sort_ex(subsectMem,ID_SlpRaw));
    Agg_data(i,ID_ElDiff)=sum(raw_data_sort_ex(subsectMem,ID_ElDiff));
    Agg_data(i,ID_Length)=sum(raw_data_sort_ex(subsectMem,ID_Length));
    Agg_data(i,ID_StrO)=(raw_data_sort_ex(subsectMem(1),ID_StrO));
    Agg_data(i,ID_MicroWSAre)=sum(raw_data_sort_ex(subsectMem,ID_MicroWSAre));
    Agg_data(i,ID_FldPlnWdth)=sum(raw_data_sort_ex(subsectMem,ID_FldPlnWdth));
    Agg_data(i,ID_Ad)=sum(raw_data_sort_ex(subsectMem,ID_Ad));
    Agg_data(i,ID_FX)=raw_data_sort_ex(subsectMem(1),ID_FX);
    Agg_data(i,ID_FY)=raw_data_sort_ex(subsectMem(1),ID_FY);
    Agg_data(i,ID_TX)=raw_data_sort_ex(subsectMem(end),ID_TX);
    Agg_data(i,ID_TY)=raw_data_sort_ex(subsectMem(end),ID_TY);
    Agg_data(i,ID_Wac)=mean(raw_data_sort_ex(subsectMem,ID_Wac));
    Agg_data(i,ID_Conf)=mean(raw_data_sort_ex(subsectMem,ID_Conf));
    Agg_data(i,23)=idSubSect;
    Agg_data(i,24)=length(subsectMem);
    % The dissolve field marks the members of each subsection. This can be used for plotting in Matlab or ArcGis 
    dissolve(subsectMem)=i;
    rand_col(subsectMem)=rand(1);
    end
 %%       
end 

Agg_data(i:end,:)=[];

close all % close previous figure
delete(h)

%% Plot the river network. 
disp('Plot aggregated river network')

    network_plotter_categories(Agg_data,Agg_data(:,24),Agg_data(:,ID_StrO),'HYMO aggregation River Network')

    