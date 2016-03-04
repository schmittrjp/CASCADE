function [ out_mat ] = parameterToMatrix(input_par,II,dir)
%PARAMETERTOMATRIX transfer a vector of source attributes into a matix of
%path attributes. E.g., write the source grain size ds in all reaches
%along the pathway from di to the basin outlet
% INPUT: 
% input_par: Vector of observed inputs 
% II: Multigraph matrix 
% dir: direction in which to write the output matrix: 
    % 1: column wise: parameter is reach dependent (e.g. slope)
    % 2: row wise: parameter is cascade dependent (e.g. grain size)
% OUTPUT: 
%out_mat: output parameter matrix with same structure as II 

out_mat=nan(size(input_par,1),size(input_par,1));

for ii=1:length(input_par) % loop through all observations 
%         waitbar(ii/length(input_par),h,'Writing parameter matrices')

%         d90ii(ii)=(hydraulicData(ii,2)); % sediment size of reach ii, derived from morphologic parameters
% %         dpercii=prctile(normrnd(d90ii(ii),0.6*d90ii(ii),100000,1), 0.1);
%         dpercii=d90ii(ii)/2.1; 

%         dii(ii)=(hydraulicData(ii,2));
% %         dii(ii)=2e-3;
%         if dpercii<=0; dpercii=0.5*10^-3; end 
%         %d10ii=prctile(normrnd(d50ii,0.2*d50ii,1000,1),0.1)
% %         hydraulicData(ii,3)=dpercii;

       jj=[find(isfinite(II(ii,:))==1)]; % loop through all reaches downstream of reach ii
%             reach_jj=find(raw_data_sort(:,ID_FromN)==jj); % position of reach j in the data
%             if numel(reach_jj)>1
%             reach_jj=reach_jj(1); %trough the resorting, some reaches are duplicate, use only the first value
%             end
            if isempty(jj)==0 % check if jj exists
              
                if dir==2
                out_mat(ii,jj)=input_par(ii); 
                elseif dir==1
                out_mat(ii,jj)=input_par(jj);    
                end 
                
                
%                 Dmat(ii,jj)=dii(ii);
%                 D90mat(ii,jj)=d90ii(ii);
%                 Thetamat(ii,jj)=2*d90ii(ii);
%                 Wmat(ii,jj)=raw_data_sort(reach_jj,ID_Wac);
%                 Lmat(ii,jj)=raw_data_sort(reach_jj,ID_Length);
                
            end % check if jj exists
        % loop through all reaches downstream of reach ii 

end % loop through all observations 

end

