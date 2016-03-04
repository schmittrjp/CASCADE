function h_selected=myLineCallback( LineH, EventData, LineList, Input, Output, PP4, Network, outlet_node_new , TranspC_mat2)

% select specific edges from a river network plot 

global h_selected cascadeResults
% disp(LineH);                    % The handle
% disp(get(LineH, 'YData'));      % The Y-data
h_selected=(find(LineList == LineH));  % Index of the active line in the list
set(LineList, 'LineWidth', 0.5);
set(LineH,    'LineWidth', 3);
set(LineH,'Color','g')
uistack(LineH, 'top');  % Set active line before all others

Path1_to_out=Network.Downstream.Path{1,1}{1,outlet_node_new(1)}; % find path from node 1 to outlet node
start_node=Path1_to_out(h_selected); % find the start node of the selected line
Path_startNode_to_out=Network.Downstream.Path{start_node,1}{1,outlet_node_new(1)}; % find the path from the start node to the outlet

inputLine=Input(start_node,Path_startNode_to_out);
outputLine=Output(start_node,Path_startNode_to_out);
tCap_line=TranspC_mat2(start_node,Path_startNode_to_out);
PPLine=PP4(start_node,Path_startNode_to_out);

plotX=1:length(Path_startNode_to_out);

figure('name', ['Sediment cascade from source ' num2str(start_node)]) 
    title(['Sediment cascade from source ' num2str(start_node)])
    plotyy(plotX,[inputLine' outputLine' tCap_line'],plotX,PPLine)
    legend('Input','Output', 'Transport Capacity', 'Percent Deposition')


end

