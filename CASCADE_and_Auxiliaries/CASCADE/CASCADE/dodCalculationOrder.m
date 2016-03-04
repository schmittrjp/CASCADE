function [ calculation_order, calculation_hierarchy] = dodCalculationOrder( AggData, Network, US_hierarchy, reservoirFromN )
% Implements downstream of Dam calculation order. The ideas is that
% cascades that begin in the river network directly downstream of a
% reservoir are routed first. This is to avoid that low order tributaries
% (that are normally routed first) that enter downstream of a resevoir take
% an over-proportional share of the transport capacity. 

% The idea is the following: 

%1: Find all nodes on the pathway from a reservoir
% to the outlet (n1). Then, find all nodes upstream or not directly
% downstream (i.e., in tributaries) of the reservoir (n2). 
%2: Determine the calculation order for n1 (co_n1)(downstream of reservoirs)
%3: Determine calculation order for all nodes (co_an)
%4: Delete nodes n1 from co_an -> this results in the calculation order for
    %the nodes n2 -> co_n2. 
%5: reassemble a new calculation order: use first co_n1 then co_n2. 

%% define some global vars 
global outlet_node_new ID_FromN

%% 1: Find all nodes n1 on the pathway from a reservoir
        n1=[];
        for res=reservoirFromN % loop through all reservoir nodes
        n1=[n1 Network.Downstream.Path{res}{outlet_node_new(1)}]; % save nodes on the pathway from the reservoirs to the outlet
        end
        
        n1=unique(n1); % Nodes might be on multiple pathways dowsntream of reservoirs. Keep each node only once.

         
%% 2: Determine the calculation order for n1 -> co_n1         
  
%2.1: Determine the hierarchy of n1  
        n1Hierarchy=US_hierarchy(n1); % hierarchy of node along the path from reservoir to the basin outlet             
     
%2.2: the calculation order of n1 is derived from sorting the hierarchy of n1
        [~,n1hierarchy]=sort(n1Hierarchy,'descend');
        co_n1=n1(n1hierarchy)';
        
%% 3: Determine calculation order for all nodes (n1 and n2)

% 3.1: Get calculation order for all nodes
        [~,co_an]=sort(US_hierarchy,'descend'); % sort the nodes according to their hierarchy in the network (e.g., node 1 is the lowest hierarchy order, the pathway node1-outlet will be calculated first). 
% 3.2: Find all nodes upstream of reservoirs (n2)
        n2=AggData(co_n1==0,ID_FromN); % all nodes upstream or not impacted by reservoirs

%% 4:  Delete nodes downstream of reservoirs  (n1) from the calculation order of all nodes
% 4.1.: Find location of n1 in co_an
        [~,n1_pos_in_co_an]=intersect(co_an,n1);
% 4.1.: Delete n1 from co_an
        co_an(n1_pos_in_co_an)=[];
       
%% 5.: Reassembling calculation order: route nodes downstream of reservoirs first, then all other nodes. 
    calculation_order=[co_n1; co_an];
    
    calculation_hierarchy(calculation_order)=1:length(calculation_order);
    
    
%% Explanation: What is the difference between calculation order and calculation hierarchy. 
% lets assume upstream of node 2421 is a dam. Hence, 2421 needs to be
% calculated first. 

% calculation hierarchy: 

% calculation_hierarchy(2421)=1. 

% Hence, node 2421 needs to be calculated first. 
% calulation order: In which order cascades need to be calculated. 
% hence:                

% calulation_order(1)=2421
       

end

