function [out] = findNDsNodes( D,S,n, Ad )
%FINDNUSNODES Finds n nodes downstream of node S along the main channel 

%%% Inputs 
%Dus: Upstream adjacency matrix (to be derived from write_adj_matrix)
%S: Start node where to begin the analysis 
%n: How many upstream nodes to be found 
%Ad: Drainage area vector, in order to find the main channel at
%bifurcations. 


nodes=zeros(n,1);

fromN=S; 

for nn=1:n % find n nodes 

    temp=D(fromN,:);
    node_pos=find(temp);
    
    
    if length(node_pos)==1 
        nodes(nn)=node_pos;
    
    elseif isempty(node_pos) % no more downstram nodes  
        nodes(nn:end)=nan; 
        break
    end 
    
    clear node_pos
    fromN=nodes(nn);
end 

out=nodes;

end

