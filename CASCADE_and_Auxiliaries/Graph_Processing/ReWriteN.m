
FN=AggData(:,ID_FromN); TN=AggData(:,ID_ToN);


newFN=(1:length(FN))'; %the IDs of the from-nodes are continous from 1 to the number of nodes. 
transferTable=[FN newFN]; % the transfer table maps these new numbers to the previous from-nodes 

% for the to-node the remapping is more difficult. 
i=0; % counter for the outletnode

for tn=unique(TN)' % loop through all to nodes 
    
    oldTN=tn ; % get the old ID of the current to node
    oldTNPos=find(TN==oldTN); % find at which position it was used (it can be at multiple positions, because a node can be a to-node for multiple reaches at confluences). 
    
    if isempty(transferTable(transferTable(:,1)==oldTN,2)) % if it was not used at all (this means that the current to  node is at the end of the network, and is therefore not found amongst the from nodes)
     i=i+1; 
     newTN(oldTNPos)=transferTable(oldTNPos,2);
     outletnode(i)=newTN(oldTNPos); % control this afterwards: If all went right there should be only one value in here. 
     else
    newTN(oldTNPos)=transferTable(transferTable(:,1)==oldTN,2);
     end 
    
end

% curNodeID=0 % current lowest ID that can be assigned to a node. 
% 
% %%
% 
% for node=1:length(FN)
%   if FN(node)>curNodeID
%        curNodeID=curNodeID+1;
%        TN(TN==FN(node))=curNodeID;
%        FN(FN==FN(node))=curNodeID;       
%        FN(node)=curNodeID;
%      
%    
%    end
%    
%    if TN(node)>curNodeID
%        
%        curNodeID=curNodeID+1;
%    
%     FN(FN==TN(node))=curNodeID;
%        TN(TN==TN(node))=curNodeID;
%        TN(node)=curNodeID; 
%        
%        
%    end
% end