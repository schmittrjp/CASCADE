clear Agg_data
% find all  headwater nodes (definition: a headwater nodde is a node that is only a from-node, but no to-node)

disp('Identifying sections')

headwaters=raw_data_sort_ex(~ismember(raw_data_sort_ex(:,ID_FromN),raw_data_sort_ex(:,ID_ToN)),ID_FromN); % very rapid approach from: http://www.mathworks.com/matlabcentral/newsreader/view_thread/300928

% find all confluences (definition: a confluence node is a node that is to-node for multiple edges.) 

[temp, ~]=histcounts(raw_data_sort_ex(:,ID_ToN),1:length(raw_data_sort_ex));
confluences=raw_data_sort_ex(temp>1,ID_FromN); clear temp

sect_start_nodes=[headwaters; confluences];

% write section membership table: which nodes belong to which section? 


for iii=1:length(sect_start_nodes)
    sectMS{iii,1}=sect_start_nodes(iii);
end 

h = waitbar(0,'Initializing waitbar...');


for ll=1:length(sect_start_nodes)
    
   waitbar(ll/length(sect_start_nodes),h,'Identifying sections')
   
   [~, sectMS{ll,2}, ]=graphshortestpath(D,sect_start_nodes(ll),outlet_node); % find  the path from each section start node to the outlet.
   temp=find(ismember(sectMS{ll,2},sect_start_nodes),2); % along that path find which nodes mark the beginning of another section
  
   if length(temp) > 1 % check if there is any downstream confluence (not the case, e.g., for reaches close to the basin outlet)
   next_confluence=temp(2); % the second value in temp is the next confluence, the first value is actually the current node
 
   temp2=sectMS{ll,2}; % get the entire path 
   temp2(next_confluence:end)=[]; % delete all nodes downstream of thenext confluence 
   sectMS{ll,2}=temp2; % save the resulting section
   end  
   
end

delete(h)
% loop through all sections and aggregate reaches that belong to the same HYMO type 

%%
h = waitbar(0,'Initializing waitbar...');

Agg_data=zeros(length(raw_data_sort_ex),22);

i=0;
for s=1:length(sectMS)

    waitbar(s/length(sectMS),h,'Aggregate subsections')
    
    sectionMem=sectMS{s,2}; % all members of a section 
    hymoClMem=HymoClassRes.Cid(sectionMem); %HYMO class of section members
    
    % detect every downstream change in HYMO class
    type_change=abs(diff(hymoClMem)); 
    temp=type_change>0;
    type_change=[1 temp]; clear temp; 
    
    % every group of nodes belonging to the same HYMO class
    % gets its own subsection
   
    subsectionID=cumsum(type_change);
   
    %%
    for ss=unique(subsectionID)
        i=i+1;
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
    Agg_data(i,23)=mean(HymoClassRes.Cid(subsectMem));
    
    % The dissolve field marks the members of each subsection. This can be used for plotting in Matlab or ArcGis 
    dissolve(subsectMem)=i;
    end
 %%       
end

Agg_data(i:end,:)=[];

close all % close previous figure
delete(h)

%% Plot the river network. 
disp('Plot aggregated river network')

    network_plotter_categories(Agg_data,Agg_data(:,23),Agg_data(:,ID_StrO),'HYMO aggregation River Network',[])
