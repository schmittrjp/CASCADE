function [h_net] = networkPlotter(raw_data, varargin)
% Plots the river network and visualizes CONTINOUS DATA. . 

% For  Can use both line width and color to represent
% different attributes. Colorcode can be clarified using a legend. 

%%% Inputs
% 1.: raw_data: Matrix containing from- and to-node infomation and attribute
% values. 
% 2.: c_att: Vector of observations of the attribute that is shown by color
% code 
% 3.: c_class: color classes to be considered
% 4.: w_att: Vector of observations of the attribute that is shown by line width
% 5.: cMap: Color map to be used 
% 6: fig_name: Name of the output figure. 
% 7.: c_att_name: Description of color attribute. If used, a legend is added to
% the figure. 
% 8.: categorial: 1 or 0; if 1, a predefined color categorie is used for each
% reach. (e.g. strahler Order, cluster ID etc.)

disp('Plotting river network')


    if isempty(varargin{5})==0 
    fig_name=varargin{5};
    figure('Name',fig_name)
    end

global lcl ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY ID_Wac ID_Q15 outlet_node% Clear temporary variables

catInd=varargin{7}; % load the switiching variable for the mode


%store the handles for each line 
hh=zeros(length(raw_data),1);

%% plot the river network 
 if isempty(varargin{1})==0 && isempty(varargin{3})==0 % Plot with linewidth and color attribute
    cAtt=varargin{1}; cClass=varargin{2}; wAtt=varargin{3};
         if size(cAtt,1)>size(cAtt,2); cAtt=cAtt'; end; % cAttClassMem must be a row vector to make th esubsequent for loop work
  
%     cmap=flipud(cmap);
% Plot a grey line if any observations has nan values 
    
   for ll=find(isnan(cAtt))
            hh(ll)=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)],...
            'color',[0.50 0.5 0.5],'linewidth',wAtt(ll)/2);

   end
 
 if catInd==0; % input data not categorial, but continous
    cMapLength=length(unique(cClass))+1;
    cMapName=[varargin{4} '(' num2str(cMapLength) ')'];
    cmap=eval(cMapName);
%     cmap=parula(length(unique(cClass))+1); % color map definition
    
  % loop through all classes 
    for c_cl=1:length(cClass)+1
      c_cl  ;
      % find all observations that have an attribute value that falls
      % within the current c class
      if c_cl==1; 
          cClassMem=find(cAtt<=cClass(c_cl)); 
      elseif c_cl<=length(cClass) && c_cl>1
           cClassMem=find(cAtt>cClass(c_cl-1) & cAtt<=cClass(c_cl));
      else
          cClassMem=find(cAtt>cClass(c_cl-1));
      end 
%       catt(catt==0)=nan; 
          
        for ll=(cClassMem) %=find(raw_data(:,ID_FromN)>0)'; % plot entire river network 
           
            if ll==1101
            debug=1; 
            end
 
            hh(ll)=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)],...
            'color',cmap(c_cl,:),'linewidth',wAtt(ll)/2);

        end
    
    end
 else
     
%        cmap=cbrewer('qual','Set1',length(unique(cClass))); % color map definition
%         cmap=parula(length(unique(cClass)));
    cMapLength=length(unique(cAtt));
    cMapName=[varargin{4} '(' num2str(cMapLength) ')'];
    cmap=eval(cMapName);  
    cmap=cmap(randperm(length(cmap)),:); % randomly permutate cmap for categorial data (hence without a natural gradient in data)
    
       
         for ll=1:length(cAtt) %=find(raw_data(:,ID_FromN)>0)'; % plot entire river network 

          if isnan(cAtt(ll))==0 && cAtt(ll)>0 % check if the current reach has a color categorie assigned
            c_cl=cAtt(ll);
            hh(ll)=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)],...
            'color',cmap(c_cl,:),'linewidth',wAtt(ll)/2);
        
          else %if not, plot it gray
            hh(ll)=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)],...
            'color',[0.3 0.3 0.3],'linewidth',wAtt(ll)/2); 
          end
          
          
        end
 end   

  
 %% no color or width attribute given: Just plot the "empty" river network. 
 
 elseif isempty(varargin{1}) && isempty(varargin{2})
     
    for ll=find(raw_data(:,ID_FromN)>0)'; % plot entire river network 
            hh(ll)=line([raw_data(ll,ID_FX) raw_data(ll,ID_TX)],[raw_data(ll,ID_FY) raw_data(ll,ID_TY)]);

    end
     
 end 
 

 
 
    
 
 
%% if required plot a legend for the color attribute (see list of input args)
if isempty(varargin{6})==0;
   hold on 
   
   %create fake lines for the legend. 
   for leggg=1:length(cClass)
         hh_leg(leggg)=line([raw_data(1,ID_FX) raw_data(1,ID_FX)],[raw_data(1,ID_FY) raw_data(1,ID_FY)],'color',cmap(leggg,:),'linewidth',2);  % for each color attribute value make a fake line with length 0. 
   end 
           if max(cClass)>10^1; legEnt=num2str(cClass','%10.2E'); else legEnt=num2str(cClass'); end; % prepare the text for the legend

   % create legend 
%    leg=legend(hh_leg,legEnt,'Location','east');
   leg=legend(hh_leg,legEnt);
    set(leg,'Position', [0.575    0.6    0.2232    0.2214]);
   set(leg,'FontSize',12);
   leg_pos=get(leg,'Position');
   text(leg_pos(1), leg_pos(2)+leg_pos(4)+0.1,varargin{6},'Units','normalized','Fontweight','bold','FontSize',16);
   
end 

scalebar('km');

   axis equal off
h_net=hh; 


set(h_net, 'ButtonDownFcn', {@LineCallbackNetwork, h_net});

   
end



