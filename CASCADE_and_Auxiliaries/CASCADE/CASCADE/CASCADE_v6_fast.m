%%%%%%%% CASCADE core model for sediment routing along multiple sediment cascades. 

%%% Changed calculataion of competition for improved efficiency

%%% Changed calculation of QSij' (competition corrected transport capacity)
%%% if there are reservoirs. Now, QSij is recalculated after the
%%% construction of reservoirs. Sources downstream of reservoirs can
%%% provide more sediment, now. 

%%% Changed calculation of reservoir impacts on tcap. TrnspC_mat1 is now
%%% adapted to consider for reservoirs before the running of the model.
%%% This enables the use of different routing schemes.

%%% Changing source supply can now be specified as maximum incision rate. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define calculation order of cascades
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist('calculationOrderType')==0; calculationOrderType='upstream'; end ; % check if calculation order is defined or not
    
    if strcmpi(calculationOrderType,'upstream') %upstream nodes routed first (standard)
    [~,calculation_order]=sort(US_hierarchy,'descend'); % sort the nodes according to their hierarchy in the network (e.g., node 1 is the lowest hierarchy order, the pathway node1-outlet will be calculated first). 
    elseif strcmpi(calculationOrderType,'downstream') %downstream nodes routed first
    [~,calculation_order]=sort(US_hierarchy,'ascend');
    elseif strcmpi(calculationOrderType,'dod'); %dod: Downstream of DAMs. Cascades Downstream of dams are routed first. Like this, these cascades can supply more sediment to simulate incision 
    
       % this calculation order makes only sense if there are reservoirs defined, and the reservoir option is activated
       if strcmpi(addReservoir,'yes') & isempty(reservoirFromN)==0
        [calculation_order, dod_calculation_hierarchy]=dodCalculationOrder( AggData, Network, US_hierarchy, reservoirFromN );
       else 
           warning('Downstream-of-dam routing scheme selected, but no reservoirs are defined or reservoir option not activated. Using standard routing scheme instead');   
                [~,calculation_order]=sort(US_hierarchy,'descend'); % sort the nodes according to their hierarchy in the network (e.g., node 1 is the lowest hierarchy order, the pathway node1-outlet will be calculated first). 
       end 
        % plot of the calculation order.    
%       h_net=networkPlotter(AggData,dod_calculation_hierarchy,1:3700, AggData(:,ID_StrO),'parula','Calculation hierarchy',[],0);
    
    elseif sum(strcmpi(calculationOrderType,{'downstream', 'downstream', 'dod'}))==0 %downstream nodes routed first
               warning('calculationOrderType not valid. Using standard (uptream cascades first) routing scheme instead');   
                [~,calculation_order]=sort(US_hierarchy,'descend'); % sort the nodes according to their hierarchy in the network (e.g., node 1 is the lowest hierarchy order, the pathway node1-outlet will be calculated first). 

    end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate REFERENCE TRANSPORT CAPACITY (i.e., for the d50). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decide of the source grain sizes are to be weighted by a certain factor
weightMat=Dmat./Dmat; 
    for rr=1:length(TranspC_mat1)
        c_in=find((TranspC_mat1(:,rr))>0); % incoming cascades into reach rr
        if isempty(c_in)==0
        d50(rr)=weightedMedian(Dmat(c_in,rr),weightMat(c_in,rr));
        end
    end
    
% Calculate the transport capacity for the local d50   
calculateQsd50

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate REFERENCE TRANSPORT CAPACITY (i.e., for the d50). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TranspC_mat1=Qsij; %; TranspC_mat1(TranspC_mat1==0)=nan; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize storage structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Tmat=nan(size(TranspC_mat1));
   TranspC_mat2=nan(size(Tmat));
   TranspC_mat3=TranspC_mat1;
   TranspC_mat2_noRes=nan(size(TranspC_mat1)); %% new in this script: Reference transport capacity without reservoirs
   SourceSupplyLimit=nan(size(Tmat));
   F1=nan(size(Tmat));
   F1_noRes=nan(size(Tmat)); %% new in this script: competition factor without reservoirs
   Diffmat=nan(size(Tmat));
   Input=nan(size(Tmat));
   Output=nan(size(Tmat)); 
   Distmat_cum=nan(size(Tmat)); 
   ncon_reaches=nan(size(Tmat)); 
   Discon_for_mbal=zeros(size(Tmat,1),3);
   Discon_for_local=zeros(size(Tmat,1),3);
   free_tcap=zeros(1,size(Tmat,1));
   jj=0;
   
   Input_0=diag(TranspC_mat1); % Input in the first reach of all pathways
   Input_1=zeros(size(Input_0));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correct for the presence of reservoirs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% If there are reservoirs, there need to be 2 adaptations. 
% 1: Cascades that begin upstream of the reservoir will be not present any
% more downstream of the reservoir. Hence, they will also not participate
% in competition. If cascade c_i begins upstream of the reservoir, all
% values of c_1 donwnstream of the reservoir will be deleted from the
% Transport capacity matrix
 
   if exist('addReservoir','var')==0; addReservoir='No'; end
   
        switch addReservoir  % check if higher level script asks for adding reservoirs 
                      
            case 'Yes' % if yes: a reservoir is represented as a reach with 0 transport
                % capacity, and the current cascade will not be present in
                % the river network downstream of the reservoir 
                
                for ii=calculation_order' 
                ds_path_nodes=Network.Downstream.Path{ii}{outlet_node_new(1)}; % find all downstream nodes in the proper order along the path       
 
                    %%% check if there are any nodes on the current pathway downstream of a reservoir    
                    if sum(intersect(ds_path_nodes,reservoirFromN))>0 

                        % if yes: reservoirDSnodes defines which nodes are actually downstream of a reservoir.  
                        for rdn=intersect(ds_path_nodes,reservoirFromN)
                        reservoirDSii=find(ds_path_nodes==rdn); % where is the location of the reservoir  
                        TranspC_mat1(ii,ds_path_nodes(reservoirDSii:end))=nan; % delete transport capacity downstream of the reservoir 
                        end

                    end; clear ii
                
                end
%                 TranspC_mat2(:,reservoirFromN)=nan;                                
                  
        end    
           
% 2.: Sources dowsntream of reservoirs migh produce more sediment. Check if this increase is limited from a higher-order script and display the selected value.      
%Source supply limit states how much more sediment a source can produce per year in case there is an upstream reservoir.       
               
      if exist('maxTcapIncrease','var')==0; 
         maxTcapIncrease=0; disp(['maxSupplyIncrease Variable not defined.  Source supply limit = 1']); 
      else 
         disp(['Source supply limit = ' num2str((maxTcapIncrease+1)*100) ' %'])
      end 
    
      if exist('maxIncision','var')==1; 
         disp(['Maximum incision Rate is ' num2str(maxIncision) ' m/yr']); 
      end 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the CASCADE model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=waitbar(0,'Running CASCADE');

for ii=(calculation_order)' %loop through all sources/cascades 
 
    jj=jj+1; % running counter for cascades

    waitbar(jj/length(calculation_order),h);
    
    clear ds_path_nodes_tcap; rd_discon=0;  flag=0; % reset some variables 

 %%%% Check if source ii is topologically connected at all to the basin outlet
    if isempty(Network.Downstream.Path{ii}{outlet_node_new(1)})==1
        Tmat_cum(ii,:)=0;
        Pmat_cum(ii,:)=0;
        Diffmat(ii,:)=0;

    %%%% if source ii is connected:  extract the pathway to the outlet  
        else %   
        ds_path_nodes=Network.Downstream.Path{ii}{outlet_node_new(1)}; % find all downstream nodes in the proper order along the path       
%         if ii==3115
%             asd=0; 
%         end 

    %%%% Identify nodes downstream of ii that are not connected because dii canot be transported or because there is areservoir, and update pathway accordingly. 
    % this steo is important for the computational effort. Disconnected
    % nodes can be excluded a-priori from the calculations.  
                
            voidnotes=find(isfinite(TranspC_mat1(ii,ds_path_nodes))==0 | TranspC_mat1(ii,ds_path_nodes)==0,2,'first'); % find which of the downstream nodes are not sedimentologically connected
            voidnotes2=ds_path_nodes(find(isfinite(TranspC_mat1(ii,ds_path_nodes))==0 | TranspC_mat1(ii,ds_path_nodes)==0)); % absolute ID of the disconnected downstream nodes        
            
            
            if length(voidnotes)>=2  
              ds_path_nodes(voidnotes(2):end)=[]; % clear not connected notes from the pathway. IMPORTANT: the first void note is the last node to which sediment is deliverd. So DO NOT clear it....
              % prepare the storage for disconnectivity: 
              Discon_for_local(ii,1)=ii; % 1: which pathway was interrupted
              Discon_for_local(ii,2)= ds_path_nodes(end); % 2: n which edge the pathway was interrupted
              Discon_for_local(ii,3)=voidnotes(1); % 3: how long was the pathway before it was interrupted  
            end
            
            % The cascade does not participate in competition downstream of
            % the interruption. 
            TranspC_mat1(ii,voidnotes2)=nan; 
            TranspC_mat1_noRes(ii,voidnotes2)=nan; 
            
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate competition for the current cascade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate competition and competition corrected transport capacity for the current cascade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                switch scenario        
                 case 'Scenario 1' % proportional to local transport capacity 

                   F1(ii,ds_path_nodes)=TranspC_mat1(ii,ds_path_nodes)./nansum(TranspC_mat1(:,ds_path_nodes));                 
                   TranspC_mat2(ii,ds_path_nodes)=Qsd50(ds_path_nodes).*F1(ii,ds_path_nodes); 
                   
                 case 'Scenario 2' % proportional to the initial input into a reach 

                        for kk=1:length((Qsd50)) %loop through all reaches 
                            incomming_reaches=find(TranspC_mat1(:,kk)>0 & isfinite(TranspC_mat1(:,kk))); % First, find all pathways that are connected to any reach kk
                            F2(incomming_reaches,kk)= Input_0(incomming_reaches)./nansum(Input_0(incomming_reaches)); 
                        end 
                    TranspC_mat2=bsxfun(@times,Qsd50,F2); %F2.*Qsd50;  

                 case 'Scenario 3' % proportional to the source input 
                     
                    F1(ii,ds_path_nodes)=TranspC_mat1(ii,ds_path_nodes)./nansum(TranspC_mat1(:,ds_path_nodes));

                    % Define if the competition corrected or the reference transport capacity at the source is used in Scenario 3.   
                
%                   Input_1(ii)=Input_0(ii)*F1(ii,ds_path_nodes(1)); % calculate competition corrected tranport capacity in the first node of the cascade: -> initial sediment supply. 
                    
                    Input_1=Input_0; % use reference transport capacity
                    
                    TranspC_mat2(ii,ds_path_nodes)=Input_1(ii).*F1(ii,ds_path_nodes);   % actual transport capacity (with Reservoirs)
                    
                        
                end
%                 
           switch addReservoir  % check if higher level script asks for adding reservoirs 
           
           case 'Yes'  % Check how much more sediment sources supply after the construction of reservoirs and if that is within the given limit. 
               % in general a limit should be given, 
            
            % this does only work if a pre-disturbance transport capacity is defined.    
            if exist('preDisturbance','var')==0; error('Reservoirs defined but no pre-disturbance state available'); break; end;     
            

                    % calculate source supply limit for current cascade
                    SourceSupplyLimit(ii,ds_path_nodes)=preDisturbance.TranspC_mat2(ii,ds_path_nodes).*(1 + maxTcapIncrease); 
                    correct_sources=find(TranspC_mat2(ii,ds_path_nodes)>SourceSupplyLimit(ii,ds_path_nodes));
                    
                        % this option is activated if a maximum incision
                        % rate at the source reach in terms of m/yr was
                        % defined. 
                        if exist('maxIncision')         
                            VIncMax(ii)=Lmat(ii,ii)*Wmat(ii,ii)*maxIncision; % maximum volumetric rate of incision (m3/yr)
                            QIncMax(ii)=VIncMax(ii)*2600;% maximum load resulting from that incision. 
                            
                            if QIncMax(ii)<SourceSupplyLimit(ii,ii) % check if the incision limit is below the transport capacity limit. 
                            SourceSupplyLimit(ii,ii)=QIncMax(ii); % change the initial supply to be a function of incision rate. 
                            end
                        end 
                    
                    
                    % correct sources where supply limit is exceeded
                    if length(TranspC_mat2(ii,ds_path_nodes)>SourceSupplyLimit(ii,ds_path_nodes))>1
                    TranspC_mat2(ii,ds_path_nodes(correct_sources))=SourceSupplyLimit(ii,ds_path_nodes(correct_sources));        
                    end  
                    
%                     SourceSupplyLimit(ii)=preDisturbance.TranspC_mat2(ii,ii).*(1 + maxSupplyIncrease); 
%                     if  TranspC_mat2(ii,ii)>SourceSupplyLimit(ii)
%                         TranspC_mat2(ii,ii)= SourceSupplyLimit(ii);                         
%                     end
                     

           end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Routing of the current sediment cascade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tcap=TranspC_mat2(ii,ds_path_nodes); % transport capacity along the cascade                    
            % prepare storage 
                input=zeros(size(tcap));
                output=zeros(size(tcap));
                delta=zeros(size(tcap));
           
            %loop through all edges on cascade ii
            for rd=1:length(tcap) 
    
                if rd==1 %first upstream edge, no inputs, only outputs 
                    
                    input(rd)=0; 
                    output(rd)=tcap(rd);
                    delta(rd)=input(rd)-output(rd);
                     
                else % if second or  more reach along the path, calculate the sediment balance 
                      
                    input(rd)=output(rd-1); % the input equals the output from the previous reach 
                    
                    if isnan(input(rd)) || input(rd)==0 % checkk if there is actually any input coming in (if not: stop)
                       
                    flag=1; % no input, e.g. because of reservoir: break 
                        
                    elseif input(rd)<=tcap(rd)  % input smaller than transport capacity 
                    output(rd)=input(rd);   % -> input=output, all sediment is routed through the reach 
                    delta(rd)=input(rd)-tcap(rd);

                    else                    % input is larger than tranport capacity 
                    output(rd)=tcap(rd);    % -> output equals the transport capacity     
                    delta(rd)=input(rd)-output(rd);
                    end    

                    % if a threshold value of deposition is exceeded, interrupt the sediment transport 
                       if output(rd)/output(1) < transport_treshold && isnan(output(rd))==0; 
 
                            output(rd)=nan; 
                            
                            % update transport capacity matrices if a cascade
                            % is interrupted. This cascade will not partake any more
                            % in competition. 
                     
                            TranspC_mat1(ii,ds_path_nodes(rd))=nan;
                            TranspC_mat1_noRes(ii,ds_path_nodes(rd))=nan;
                            TranspC_mat2(ii,[ds_path_nodes(rd)])=nan;

                            ds_path_nodes_tcap=ds_path_nodes(rd:end); % nodes that are disconnected from the current cascade because of competition
                            
                            rd_discon=rd_discon+1; % counter for disconnectivity 
                            
                                % prepare another storage for disconnectivity (due to deposition rather than local morphologic controls): 
                                if rd_discon==1 
                                Discon_for_mbal(ii,1)=ii; % store that the sediment path was interrupted because of mass transport rather than time limitation
                                Discon_for_mbal(ii,2)=ds_path_nodes(rd); % store at which node that path was interrupted
                                Discon_for_mbal(ii,3)=rd; % store after how many nodes the path was interrupted
                                end
                           
                            flag=1; % flag for stopping calculations because this pathway was broken due to competition
                            break 
                       end % if output(rd)/output(1) < transport_treshold && isnan(output(rd))==0;  
                     
               end % if rd==1 %first upstream edge, no inputs, only outputs 
                
                if flag==1 % stopping calculations along that pathway because broken due to competition
                flag=0; %reset the flag 
                break 
                end 
            end % for rd=1:length(tcap)  

            
        Input(ii,ds_path_nodes)=input;
        Input_2(ii,ds_path_nodes(2:end))=input(2:end);
        Output(ii,ds_path_nodes)=output;   
        Diffmat(ii,ds_path_nodes)=delta;
        Distmat_cum(ii,ds_path_nodes)=(II(ii,ds_path_nodes));

        if exist('ds_path_nodes_tcap') % check if there are any nodes where connectivity was interrupted due to transport capacity limitation
             % if yes, delete the nodes on this pathways downstream from the interruption from the transport capacity and the travel time matrix. 
          TranspC_mat1(ii,[ds_path_nodes_tcap voidnotes2])=nan;
          TranspC_mat1_noRes(ii,[ds_path_nodes_tcap voidnotes2])=nan;
          TranspC_mat2(ii,[ds_path_nodes_tcap voidnotes2])=nan;
   
        end   
%         
%             if ii==1726
%             disp('hello')
%             end
%         
    end

%%% TranspC_mat3 can be used to consider only cascades that have allready been routed in the competition calculations   
%        TranspC_mat3(isnan(Output))=nan; 
%        TranspC_mat3(~isnan(Output))=TranspC_mat1(~isnan(Output));
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up the disconnectivity storage matrices. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Discon_for_mbal(Discon_for_mbal(:,2)==outlet_node_new(1),:)=[]; % dont consider the natural end of sediment paths at the basin outlet
   Discon_for_local(Discon_for_local(:,2)==outlet_node_new(1),:)=[]; % dont consider the natural end of sediment paths at the basin outlet
  
   Discon_for_mbal(Discon_for_mbal(:,2)==0,:)=[];
   Discon_for_local(Discon_for_local(:,2)==0,:)=[];

close(h)    

%% calculate the travel times in each reach based on the flux    

Flux_Area_mat=Thetamat.*Wmat; % Cross sectional area of sediment flux in each reach [m2]. 
Tmat=Lmat./(Output./2600./Flux_Area_mat); % Calculate travel time of sediment in the reaches by dividing the flux by the flux (=Output) by the characteristic crossectional area of flux in each reach. 
Vmat=(Output./2600./Flux_Area_mat);

% Create the cumulative travel time matrix. It saves the cumulative travel
% time by getting the node order of each path, and than using this order to
% sum the travel times. 

for hh=(calculation_order)'   % loop through all from nodes
    ds_path_nodes=Network.Downstream.Path{hh}{1,outlet_node_new(1)};  % get the sedimet path from the from-node to the outlet 
    Tmat_cum(hh,ds_path_nodes)=cumsum(Tmat(hh,ds_path_nodes)); % sum up the travel times along the previously defined path. Tmat contains "inf" where  the sediment path is interupted. This is considered in Tmat_cum. 
end   

%% Deposition Ratio
PP4=(bsxfun(@rdivide, Output', max(Output')))'; % create a new matrix for sediment contribution 
PP4=1-PP4; % how much (cummulative percentage) of the original input in a path is deposited downstream


