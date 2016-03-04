function [HymoClassRes] = hymo_classification( raw_data, som_vars,clust_par, plot )
% Performs a fast classification based on selected hydromorphorphologic
% parameters using self-organizing maps and fuzzy clustering (see Schmitt et al 2014). 

%%% Inputs: 
% raw_data: Matrix will all HYMO parameters for all nodes. 
% som_vars: Which of these input parameters are used for the classification
% (in which columns of the input matrix these parameters are stored)
% plot: activates plotting of the SOM and of the classified river network.
% sort_var: The clusters can be sorted (e.g. such that the cluser
% representing nodes with high drainage area also has a high number). This
% variable gives which of the som_vars is used for this purpose. 

%%% Outputs: 
% HymoClassRes: Struct with most relevant results of the classification 
% HymoClassRes.Bmu: BMU membership of each node 
% HymoClassRes.Cid: Dominant cluster for each node 
% HymoClassRes.U: Fuzzy membership for eahc node

global ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY ID_Wac ID_Q15% Clear temporary variables


%% run the SOM algorithm
SOM_data=raw_data(raw_data(:,ID_FromN)>0,som_vars); %prepare input data  
               
Vars={'Slope','Flodplain W', 'Q 1.5', 'WAC'}; %prepare vector with variable names 

sD=som_data_struct(SOM_data); %make SOM input struct 

sD=som_normalize(sD,clust_par.normalization); %normalize input data (very important!!)

sM=som_make(sD,'algorithm','batch','training','long'); % run SOM algorithm 

BMU_membership=som_bmus(sM,sD); % finds to which BMU a node belongs


 % h=som_show(sM);
%% Fuzzy Clustering 

% find the optimal cluster number with the subclustering approach
[center_sc ~]=subclust(sM.codebook,clust_par.epsilon);
n_cluster=size(center_sc,1);

% print([num2str(n_cluster) 'HYMO classes identified']) 

% run the fuzzy clustering
[center_fcm U obj_fcn]=fcm(sM.codebook,n_cluster,[nan, nan, nan, 0]);

% resort the cluster centres 
% by which input variable to resort the clusters. 
[~,sort_ind]=sort(center_fcm(:,clust_par.sort_var));

U=U(sort_ind,:);

[C Fuzz_ind]=max(U); %Fuzz_ind identifies to which the cluster each BMU has the highest degree of belonging

% now, the BMUs are clustered. Now define for each node to which cluster it
% belongs. 

cid=zeros(size(SOM_data,1),1);
u=zeros(size(SOM_data,1),size(U,1));

for nn=1:size(SOM_data,1)
   
   bmu_nn=BMU_membership(nn); % find to which BMU the current node belongs (bmu_nn)
   cid(nn)=Fuzz_ind(bmu_nn); % find the dominant cluster for bmu_nn and assign it to the current node
   u(nn,:)=U(:,bmu_nn); % find the fuzzy signature for bmu_nn
   
end

clear U

consideredNodes=find(raw_data(:,ID_FromN)>0); % this is required in case some nodes, e.g., with small Strahler Order are excluded. 
% prepare outputs struct
HymoClassRes.Bmu(consideredNodes)=BMU_membership; 
HymoClassRes.Cid(consideredNodes)=cid; 
HymoClassRes.U(1:size(u,2),consideredNodes)=u'; 


% Plot results 
if nargin==4
    % Plot the SOMs
    figure('Name','Self Organizing Maps')
    h=som_show(sM);
    % Plot the classified river network
    network_plotter_categories(raw_data(consideredNodes,:),cid,raw_data(consideredNodes,ID_StrO),'Classified River Network','Hymo Class')
end
   
end