function resLocation=myLineCallback( LineH, EventData, LineList, VArs)

% select specific edges from a river network plot 

global CascadeResults outlet_node_new resLocation

v2struct(Vars)

resLocation= [resLocation (find(LineList == LineH))];  % Index of the selected edge in the list

hold on 
scatter(AggData(resLocation,ID_FX),AggData(resLocation,ID_FY),200,'>','k','filled')

% Define some options: 
          scenario='Scenario 3'; 
          Competition=3; % 1: no competition, transport capacity is infinite (all inputs transported downstream), 2: no competition, but local transport capacity limits, 3: competition between pathways, and local capacity limits 
          Tcap_redistribution=0; % 1.: free trasnport capacity can be redistributed between pathways 
          transport_treshold=1e-2; % after how much percent of the initial input a pathway is interrupted 
          addReservoir='Yes';

reservoirFromN=resLocation;
          
CASCADE_v5_fast



end

