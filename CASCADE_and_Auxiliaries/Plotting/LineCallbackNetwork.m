function h_selected=myLineCallback( LineH, EventData, LineList)

% select specific edges from a river network plot 

global CascadeResults outlet_node_new
h_selected=(find(LineList == LineH));  % Index of the selected edge in the list
% set(LineList, 'LineWidth', 0.5);
% set(LineH,    'LineWidth', 3);
% set(LineH,'Color','g')

%%%% use a scatter to mark the selected edge
startXY=[get(LineH,'Xdata'); get(LineH,'Ydata')]; % get coordinates 
hold on
scatColor=rand(1,3);
scatter(startXY(1,1),startXY(2,1),100,scatColor,'filled','LineWidth',1.5,'MarkerEdgeColor','k')

%%%% plot the pathway from that edge until it is interrupted
Path_to_out=CascadeResults.Network.Downstream.Path{h_selected}{outlet_node_new(1)}; % find the path to the outlet
interruptionNode=find(isnan(CascadeResults.Output(h_selected,Path_to_out)),1,'first')+1; % find the first node where the path is interrupted
if isempty(interruptionNode) interruptionNode=length(Path_to_out); end % if the pathway is not interrupted plot tranport to the basin outlet

interruptionNodeID=(Path_to_out(interruptionNode)); % convert the node number in the path into an absolute ID

hInterruptionNode=LineList(interruptionNodeID); % handle of interruption node
startXY=[get(hInterruptionNode,'Xdata'); get(hInterruptionNode,'Ydata')]; % get coordinates of interruption node
scatter(startXY(1,1),startXY(2,1),100,scatColor,'d','filled','LineWidth',1.5,'MarkerEdgeColor','k') % plot interruption node

%%% plot the pathway between start node and interruption node
pathX=zeros(interruptionNode,2); %storage for path x coordinates
pathY=zeros(interruptionNode,2); % storage for path y coordinates

    for pw=1:interruptionNode % loop through all nodes on the pathway
        hpw=LineList(Path_to_out(pw)); % find handle of next node on the pathway
        startXY=[get(hpw,'Xdata'); get(hpw,'Ydata')];
        pathX(pw,:)=startXY(1,:);
        pathY(pw,:)=startXY(2,:);  
    end
   hold on
   hpth=line(pathX',pathY','Color',scatColor,'LineWidth',10); % plot the pathway between the start node and the interuption node
   uistack(hpth, 'bottom');  % send that line to the bottom

%%%% plot input, output, and delivery ratio in a separate plot   
CascadeResults.Network.Downstream.Path{h_selected}{outlet_node_new(1)};
    inputLine=CascadeResults.Input(h_selected,Path_to_out(1:interruptionNode));
    outputLine=CascadeResults.Output(h_selected,Path_to_out(1:interruptionNode));
    tCap_line=CascadeResults.TransportCapacity(h_selected,Path_to_out(1:interruptionNode));
    PPLine=CascadeResults.Delivery(h_selected,Path_to_out(1:interruptionNode));
    PPLine(end)=1;

    
    plotX=CascadeResults.Network.Downstream.Distance{h_selected}(Path_to_out(1:interruptionNode))./10^3; % x coordinates for the plot 

    figure('name', ['Sediment cascade from source ' num2str(h_selected)]); % create a figure 
        [Ax, h1, h2]=plotyy(plotX,[inputLine' outputLine' tCap_line'],plotX,PPLine); 
        set(Ax(2),'Ylim',[0 1])
        % make axis labels.
        ylabel(Ax(1),'Input, Output, Trasnport Capacity [kg/yr]')
        ylabel(Ax(2),'Deposition Ratio')
        xlabel(Ax(1), 'Downstream Distance [km]')
        % set line attributes
        set(h1(1),'Linestyle',':','linewidth',1.5)
        set(h1(2),'Linestyle','--','linewidth',1.5)
        set(h1(3),'Linestyle','-.','linewidth',1.5)
        set(h2,'linewidth',2) % thicker line for deposition
        % make legend 
        legend({'Input','Output', 'Transport Capacity', 'Deposition Ratio [ ]'},'Fontsize',12);
        %make title
        title(['Sediment cascade from source ' num2str(h_selected)],'color',scatColor); % create a title 
 
   %%% Add source grain size go that plot 
       
       AxPos=get(gca,'Position');
       ax2=axes('Position',AxPos);
       dh=CascadeResults.Dmat(h_selected,h_selected);  
       scatter(ax2,0,0,sqrt(dh*50000),[0.3 0.3 0.3],'filled','MarkerEdgeColor','k','LineWidth',2)
       
       set(ax2,'xlim',[0 1],'ylim',[0 1],'visible','off')
       text(0, -0.05,['Source Grain size: \newline' num2str(dh*1000) ' [mm]'],'fontweight','bold')
end

