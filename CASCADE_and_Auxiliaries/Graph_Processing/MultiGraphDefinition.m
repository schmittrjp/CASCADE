%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script for multi-graph preparation is identical with the section from "Master Prprocessing" 
%   This script is required because the multi graph is not stored if the
%   workspace is saved.
%   RUN THIS SCRIPT if IS CASCADE IS RUN FROM PREVIOUS RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now that the number of raches is reduced it is feasible to identify
% the path between all reaches. 

% Write Adj. Matrix for aggregated network 
[DAgg,~]=write_adj_matrix(AggData(:,ID_FromN),AggData(:,ID_ToN),AggData(:,ID_Length));
[DUsAgg,~]=write_adj_matrix(AggData(:,ID_ToN),AggData(:,ID_FromN),AggData(:,ID_Length));

S=AggData(:,ID_FromN);

%vectorized approach from:  http://www.alecjacobson.com/weblog/?p=3868#respond
% this approach is elegant and finds all downstream nodes for all
% from-nodes.  
 
[Network.Downstream.Distance, Network.Downstream.Path, Network.Downstream.Predecessors] = ...
    (arrayfun(@(fromnode) graphshortestpath(DAgg,fromnode),S,'UniformOutput',false));

% Transfer downstream path from each each node into a matrix representation

II=cell2mat(Network.Downstream.Distance);
II(isfinite(II)==0)=nan;   

% find the number of downstream nodes, and the outlet-node (alldsnodes) for
% each node. The latter is required in case there are multiple outlets. 
[ndsnodes, alldsnodes]=arrayfun(@(fromnode) max(cellfun('length',Network.Downstream.Path{fromnode})),S,'UniformOutput',false);
US_hierarchy=cell2mat(ndsnodes);

% outletnode=max(unique(cell2mat(alldsnodes)));
