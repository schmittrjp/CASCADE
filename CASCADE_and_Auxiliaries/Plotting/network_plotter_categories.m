function [h_net] = network_plotter_categories(raw_data, varargin)
% Plots the river network. Can use both line width and color to represent
% different attributes. Colorcode can be clarified using a legend. 

%%% Inputs
% raw_data: Matrix containing from- and to-node infomation and attribute
% values. 
% c_att: Vector of observations of the attribute that is shown by color
% code 
% w_att: Vector of observations of the attribute that is shown by line width
% fig_name: Name of the output figure. 
% c_att_name: Description of color attribute. If used, a legend is added to
% the figure. 

if isempty(varargin{3})==0 
    fig_name=varargin{3};
    figure('Name',fig_name)
end 


global ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY ID_Wac ID_Q15 outlet_node% Clear temporary variables


%% plot the river network 
 if isempty(varargin{1})==0 && isempty(varargin{2})==0 % Plot with linewidth and color attribute
     c_att=varargin{1}; w_att=varargin{2};
     
     cmap=parula(length(unique(c_att)));

    for ll=find(raw_data(:,ID_FromN)>0)'; % plot entire river network 

        hh=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)],...
            'color',cmap(c_att(ll),:),'linewidth',w_att(ll)/2);

    end
    
 elseif isempty(varargin{1}) && isempty(varargin{2})
     
    for ll=find(raw_data(:,ID_FromN)>0)'; % plot entire river network 

        hh=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)]);

    end
     
 end 
    
 
 
%% if required plot a legend for the color attribute (see list of input args)
if isempty(varargin{4})==0;
   hold on 
   
   %create fake lines for the legend. 
   for leggg=1:length(unique(c_att))
         hh_leg(leggg)=line([raw_data(1,ID_FX) raw_data(1,ID_FX)],[raw_data(1,ID_FY) raw_data(1,ID_FY)],'color',cmap(leggg,:));  % for each color attribute value make a fake line with length 0. 
   end 
   
   % create legend 
   leg=legend(hh_leg,{num2str(unique(c_att))});
   leg_pos=get(leg,'Position');
   text(leg_pos(1), leg_pos(2)+leg_pos(4)+0.1,varargin{4},'Units','normalized','Fontweight','bold')
   
end 

   axis equal square off



end



