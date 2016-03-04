function [AggData dissolve rand_col] = SectWClassAgg(networkData,D,DUs,segPar,SegSet)
%SectWClassAgg  Implementation of a section wise classification and
%aggregation approach. A section is defined as group of reaches between two
%confluences. Reaches within each section are classified based on their
%HYMO characteristics. Subsequqent reaches of the same class are aggregated
%into morphologicall homogeneous sub-sections. 

%   Detailed explanation goes here

global ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY ID_SubWS ID_Q15 ID_Wac ID_Conf outlet_node

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% definition of sections 

    % find all  headwater nodes (definition: a headwater nodde is a node that is only a from-node, but no to-node)
    headwaters=networkData(~ismember(networkData(:,ID_FromN),networkData(:,ID_ToN)),ID_FromN); % very rapid approach from: http://www.mathworks.com/matlabcentral/newsreader/view_thread/300928

    % find all confluences (definition: a confluence node is a node that is to-node for multiple edges.) 
    [temp, ~]=histcounts(networkData(:,ID_ToN),1:length(networkData));
    confluences=networkData(temp>1,ID_FromN); clear temp

    % Section start nodes are the combination of headwater and confluence nodes 
    sectStartNodes=[headwaters; confluences];
    sectStartNodes=sort(sectStartNodes); % resort the full table of 

    % write section membership table. the section membership table is a cell array that finally identifies which nodes belong to each section.  
    for iii=1:length(sectStartNodes) % write all section start nodes. 
        sectMS{iii,1}=sectStartNodes(iii);
    end 

    h = waitbar(0,'Initializing waitbar...'); % make a waitbar for writing the section membership table 

    % use the Dijkstra algorithm to find all downstream nodes of a section
    % start node. Identify which of these nodes belong to the section.  
    for ll=1:length(sectStartNodes)
   
        waitbar(ll/length(sectStartNodes),h,'Identifying sections')

       [~, sectMS{ll,2}, ]=graphshortestpath(D,sectStartNodes(ll),outlet_node); % find  the path from each section start node to the outlet.
       temp=find(ismember(sectMS{ll,2},sectStartNodes),2); % along that path find which nodes mark the beginning of another section

       if length(temp) > 1 % check if there is any downstream confluence (not the case, e.g., for reaches close to the basin outlet)
       next_confluence=temp(2); % the second value in temp is the next confluence, the first value is actually the current section start node itself node

       temp2=sectMS{ll,2}; % get the entire path 
       temp2(next_confluence:end)=[]; % delete all nodes downstream of the next confluence 
       sectMS{ll,2}=temp2; % save the resulting section
       end  
   
    end
    
    % in rare cases there can be an empty section
    emptySections=cellfun('isempty',sectMS);
    sectMS(emptySections(:,2),:)=[];
    
    
    delete(h) % close the waitbar

%% Identification of subsections, and segmentation. 

%     % Define also connectivity in the upstream direction. This is required to
%     % smooth the data using a moving average. 
%     [DUs,~]=write_adj_matrix( raw_data(:,ID_ToN),raw_data(:,ID_FromN),raw_data(:,ID_Length));

%%

    disp('Aggregating reaches based on Hymo class')

    AggData=zeros(length(networkData),22); %AggData is the main output. It is a matrix where each segment has an individual entry with all required parameters. 
    dissolve=zeros(length(networkData),1);
    rand_col=zeros(length(networkData),1);

    h = waitbar(0); % waitbar for the segmentatio
    i=0;
    idSubSect=0; % counter for the subsections. 
    global s 
    
    for s=1:length(sectMS) % loop through all sections s

    %% Derivation of HYMO parameters for all members of section s      
    
     waitbar(s/length(sectMS),h,'Segmenting')
    
     % get all members of section s
     sectionMem=sectMS{s,2}; 
     
     if sectMS{s}==49219
         asd=12;
     end 
     
     % The measurements of HYMO parameters can be noisy, if reach length is
     % short. Therefor it is best to apply an moving average filter to filter
     % measurements. Using a filter with a filter window of length
     % n=movavgWindow, means loosing observations for n-1 nodes. Therefor, the
     % section need to be increased, to include n-1 nodes upstream of the
     % section's first, and n-1 sections downstream the section's last node. 
     
     %% Moving average to smooth observations of HYMO parameters. 
     
     % define the length of the moving average window    
     movavgWindow=floor(length(sectMS{s,2})*SegSet.resolution); % Resolution defines the length of the moving average window in function of the section length. High Resolutions means small moving average length
     if movavgWindow==0; movavgWindow=1; end; 
     
     if s==1372
         asd=1; 
     end
     % find n-1 upstream nodes    
     sectionUsNodes=findNUsNodes(DUs,sectionMem(1),movavgWindow-1,networkData(:,ID_Ad));
     
     % how to handle a situation where there are not enough upstream nodes
     % to calculate the moving average. This can happen when thee distance
     % between the start of a new section (confluence) and the headwaters
     % is small. The check is required because otherwise valid confluences
     % are deleted from the network. % There are 2 conditions to check: 
    
     if sum(isnan(sectionUsNodes)==0)<=movavgWindow-1 ...  % if there are less nodes upstream, reduce the length of the MA window
             && sum(isnan(sectionUsNodes)==0)>0; %% if the upstream section is only nan, then the entire section will be delected anyway, below. 
        movavgWindow=sum(isnan(sectionUsNodes)==0);
        sectionUsNodes=findNUsNodes(DUs,sectionMem(1),movavgWindow-1,networkData(:,ID_Ad)); % correct the length of the upstream section

     end 

     % find n-1 downstream nodes    
     sectionDsNodes=findNDsNodes(D,sectionMem(end),movavgWindow-1,networkData(:,ID_Ad));

     % Enlarge the section by these nodes. 
     sectionMem=[fliplr(sectionUsNodes') sectionMem sectionDsNodes'];
     
     % if there are no upstream nodes, delete the nans resulting from the
     % previous step
     sectionMem(isnan(sectionMem))=[];
     
     % Get the observations of HYMO parameters for all section members.
     % Which parameters are considered exactly is defined by clust_par. 
     sectionMemHymo=networkData(sectionMem,segPar); 
     
     %% identify similar nodes along section s by kmeans clustering
     
     if length(sectMS{s,2})>1 % check if the section contains multiple nodes. 
        
        % prepare input data for k-means clustering for the members of the
        % current section. 
        clustData=[sectionMemHymo];

        % normalize the input data for the clustering between 0 and 1 
        clustDataNorm=bsxfun(@rdivide, clustData, max(clustData));
        
        % calculate the relative downstream change in each parameter. 
        dsChange=abs(diff(clustDataNorm));
        
        % sum the normalized change for the parameters -> measure for how
        % strongly HYMO parameters change between reaches. This is a very important parameter to identify downstream change.  
        dsChangetot=[0; sum(abs(diff(clustDataNorm)),2)];
        
        clustDataNorm=[clustDataNorm dsChangetot];
        % check if the moving average window is not larger than the number
        % of observations
        if movavgWindow>=size(clustDataNorm,1); movavgWindow=size(clustDataNorm,1)-2; end; 
                        
        % the data are quite noisy because of the small length of
        % the reaches. Use a moving average to smooth the data 
        clustDataNormMA=tsmovavg(clustDataNorm,'s',movavgWindow,1); 
       
        % delete section members for which there is no moving average value
        % (i.e., nodes at the beginning of a section that has no upstream
        % nodes). 
%         sectionMem(isnan(clustDataNormMA(:,1)))=[];  
        sectionMem(1:movavgWindow-1)=[];
%         clustDataNormMA(isnan(clustDataNormMA(:,1)),:)=[];
        clustDataNormMA(1:movavgWindow-1,:)=[];
        
        % delete section members, that where added dowsntream after the end
        % of the original section. Length of this deletion is movavgWindow-2
        sectionMem(end-(movavgWindow-2):end)=[];
        clustDataNormMA(end-(movavgWindow-2):end,:)=[];
        
        % Define the optimal cluster number using subtractive clustering 
        
        if length(sectMS{s,2})<5 % very short reaches are all not segmented. 
        nClust=1;
        else
        [C,~]= subclust(clustDataNormMA,SegSet.epsilon); % use the subclustering algorithm. Epsilon controls helps to identify how many clusters are to be used 
        nClust=size(C,1); %number of clusters to be used
        end
                
        kCid=kmeans(clustDataNormMA,nClust);%find to which cluster center to which each node in the current section belongs. 

%         test=AbsChangeSlp(1:end-1)./AbsChangeSlp(2:end);

        % identify segments: subsequent nodes that belong to the same
        % cluster. 
        typeChange=abs(diff(kCid)); % find change in the nodes cluster.  
        temp=[1; typeChange]>0; 
        typeChange=[temp']; clear temp; 
        
    else % if there is only one node in the section.  
     typeChange=1;   
    end
    
    % every group of nodes belonging to the same HYMO class
    % gets its own segment ID 
     segmentID=cumsum(typeChange);
    
%         figure
%             subplot(2,1,1)
%             stairs(clust_data_norm_mov)
%             legend('slp','fldpln','conf','chnage')
%             subplot(2,1,2)
%             stairs(kCid);ylim([0 n_clust+1])
   %% Calculate the Hymo parameters for each segment.     
   
    %
    for ss=unique(segmentID) % loop through all segments of the current segment
        i=i+1;idSubSect=idSubSect+1; % increase the segment counter i 
        
    segmentMem=sectionMem(segmentID==ss);
    
    AggData(i,ID_arcid)=networkData(segmentMem(1),ID_arcid);
    AggData(i,ID_FromN)=networkData(segmentMem(1),ID_FromN);
    AggData(i,ID_ToN)=networkData(segmentMem(end),ID_ToN);
    AggData(i,ID_ElUs)=networkData(segmentMem(1),ID_ElUs);
    AggData(i,ID_ElUsRaw)=networkData(segmentMem(1),ID_ElUsRaw);
    AggData(i,ID_ElDs)=networkData(segmentMem(end),ID_ElDs);
    AggData(i,ID_ElDsRaw)=networkData(segmentMem(end),ID_ElDsRaw);
    AggData(i,ID_Slp)=mean(networkData(segmentMem,ID_Slp));
    AggData(i,ID_SlpRaw)=mean(networkData(segmentMem,ID_SlpRaw));
    AggData(i,ID_ElDiff)=sum(networkData(segmentMem,ID_ElDiff));
    AggData(i,ID_Length)=sum(networkData(segmentMem,ID_Length));
    AggData(i,ID_StrO)=(networkData(segmentMem(1),ID_StrO));
    AggData(i,ID_MicroWSAre)=sum(networkData(segmentMem,ID_MicroWSAre));
    AggData(i,ID_FldPlnWdth)=sum(networkData(segmentMem,ID_FldPlnWdth));
    AggData(i,ID_Ad)=sum(networkData(segmentMem,ID_Ad));
    AggData(i,ID_FX)=networkData(segmentMem(1),ID_FX);
    AggData(i,ID_FY)=networkData(segmentMem(1),ID_FY);
    AggData(i,ID_TX)=networkData(segmentMem(end),ID_TX);
    AggData(i,ID_TY)=networkData(segmentMem(end),ID_TY);
    AggData(i,ID_SubWS)=min(networkData(segmentMem,ID_SubWS));
    AggData(i,ID_Wac)=mean(networkData(segmentMem,ID_Wac));
    AggData(i,ID_Q15)=mean(networkData(segmentMem,ID_Q15));
    AggData(i,ID_Conf)=mean(networkData(segmentMem,ID_Conf));
    AggData(i,24)=idSubSect;
    AggData(i,25)=length(segmentMem);
    % The dissolve field marks the members of each subsection. This can be used for plotting in Matlab or ArcGis 
    dissolve(segmentMem)=i; % the dissolve vector can be used to dissolve to which segment an original reach belongs. 
    rand_col(segmentMem)=rand(1)*10.^(rand*5); % this random number is useful for the plotting in ArcGIS, to assign a different color to each segment, while makin sure that neighbor segents get quite different colors
     
    end
 %%       
end 

AggData(i+1:end,:)=[]; % delete empty rows in the AggData

close all % close previous figure
delete(h)

%% Plot the river network. 
% disp('Plot aggregated river network')
% 
%     network_plotter_categories(AggData,AggData(:,24),AggData(:,ID_StrO),'HYMO aggregation River Network')

    

end

