function [D,pos_tonodeID] = write_adj_matrix( fromN,toN,length )
%write_adj_matrix This function takes a two vectors of from-nodes and to-nodes transfers them into a (sparse) adjacency matrix . Each pair of from-nodes and to-nodes
%%% defines a river reach. Vectors need to be of equal size. 

%%% Input: 
% FromN: Vector of from-nodes
% toN: Vector of to-nodes


%%% Output
% D: Sparse adjacency matrix 
% pos_tonodeID: Sorting index to make sure that the position of a reach and
% its from-Node are identical (i.e., that reach 465 is on the 465th
% position in the data). 

adj_mat_ds=spalloc(max(fromN),max(fromN),max(fromN)); %(max(fromN),max(toN));
k=0


    for fromnode=1:max(fromN)
        hwb=waitbar(fromnode/max(fromN)); set(hwb,'Name','Writing Adjacency Matrix')
        fromnode;
        
        % for a given from node (rows), find to which node it is connected
        % (columns)
        %   1  2   3   4
        %1  0  23  0   0 -> from node 1 to node 2, distance is 23
        %2  0  0  34   0
        %3  0  0   0   45
        %4  0  0   0   0

        pos_fromN=find(fromN==fromnode);


       
        if isempty( pos_fromN)==0
           pos_tonodeID(fromnode,1)=pos_fromN; % store the position of each fromnode
            % store the distance 
            l=length(pos_fromN);
            adj_mat_ds(fromN(pos_fromN),toN(pos_fromN))=l;
        else
            k=k+1;
            emptynode_ID(k)=fromnode;
            pos_tonodeID(fromnode,1)=fromnode;
        end
    end
     
        
% create sparse adjacency matrix;
D=adj_mat_ds;



% S=unique(fromN);  
% [Network.Downstream.Distance, Network.Downstream.Path, Network.Downstream.Predecessors] = (arrayfun(@(fromnode) graphshortestpath(D,fromnode),S,'UniformOutput',false))

% 
% % resort data matrices to match the order of from nodes
%     raw_data_sort=raw_data(pos_tonodeID,:);
% % h = view(biograph(D,[],'ShowWeights','on'))
% 
% %% Find shortest downstream path for each node
%          S=1:max(raw_data_sort(:,ID_FromN));
% 
% % calculate connectivity matrix -> which node is connected to which
% % downstream nodes 
%   [Network.Downstream.Distance, Network.Downstream.Path, Network.Downstream.Predecessors] = (arrayfun(@(fromnode) graphshortestpath(D,fromnode),S,'UniformOutput',false));
% 
%      %vectorized approach from:  http://www.alecjacobson.com/weblog/?p=3868#respond
%      % this approach is elegant and finds all potential sinks from a given from-node: This is not actually required in a downstream river network analysis, as per definition, there is only one downstream path. 
%      
% 
%          II=cell2mat(arrayfun(@(fromnode) graphshortestpath(D,fromnode)',S,'UniformOutput',false));
%          II=II'; %II contains: in the rows the from node, and in the columns the travel distance to all nodes to which the node is connecte ddownstream
%          II(isfinite(II)==0)=nan;   
%    % read the nodes along all downstream paths 
%    [ndsnodes, alldsnodes]=arrayfun(@(fromnode) max(cellfun('length',Network.Downstream.Path{fromnode})),S,'UniformOutput',false);
%    outletnode=max(unique(cell2mat(alldsnodes)));
%    
% %% Hierarchy analysis: Find teh maximum path length UPSTREAM of each node. 
% S=1:max(raw_data_sort(:,ID_FromN));
% % calculate connectivity matrix -> which node is connected to which
% % UPSTREAM nodes  
% Dup=D'; % inverse the adjacency matrix, so the analsis is upstream. 
%   [Network.Upstream.Distance, Network.Upstream.Path, Network.Upstream.Predecessors] = (arrayfun(@(fromnode) graphshortestpath(Dup,fromnode),S,'UniformOutput',false));
%   temp=arrayfun(@(fromnode) max(cellfun('length',Network.Upstream.Path{1,fromnode})'),S,'UniformOutput',false);     
%   US_hierarchy=cell2mat(temp);  
%   
%          II=cell2mat(arrayfun(@(fromnode) graphshortestpath(D,fromnode)',S,'UniformOutput',false));
%          II=II'; %II contains: in the rows the from node, and in the columns the travel distance to all nodes to which the node is connecte ddownstream
%          II(isfinite(II)==0)=nan;   
%    % read the nodes along all downstream paths 
%    [ndsnodes, alldsnodes]=arrayfun(@(fromnode) max(cellfun('length',Network.Downstream.Path{fromnode})),S,'UniformOutput',false);
%    outletnode=max(unique(cell2mat(alldsnodes)));
%    
% [Network.Downstream.Distance, Network.Downstream.Path, Network.Downstream.Predecessors] = (arrayfun(@(fromnode) graphshortestpath(D,fromnode),S,'UniformOutput',false));
% 
% 
% save('Connectivity_data')
% 
close(hwb)
end
% 
