%%Determines a characteristic grain size for each reach, and writes these
%%char. grain sizes along all pathways. Determines also soem additionally
%%requrie dmorphologic matrices 

%OUTPUT
%1.: Dmat
%2.: Thetamat -> typical transport depth is a function of the grain size. 
%3.: Wmat
%4.: Lmat

disp('Writing parameter matrices')
h = waitbar(0,'Writing parameter matrices');

Dmat=zeros(n_obs,n_obs);
D90mat=zeros(n_obs,n_obs);
Thetamat=zeros(n_obs,n_obs);
Wmat=zeros(n_obs,n_obs);
Lmat=zeros(n_obs,n_obs);


for ii=1:size(AggData,1) % loop through all observations 
        waitbar(ii/length(AggData),h,'Writing parameter matrices')

%         d90ii(ii)=(hydraulicData(ii,2)); % sediment size of reach ii, derived from morphologic parameters
% %         dpercii=prctile(normrnd(d90ii(ii),0.6*d90ii(ii),100000,1), 0.1);
%         dpercii=d90ii(ii)/2.1; 

        dii(ii)=(hydraulicData(ii,2));
%         dii(ii)=2e-3;
        if dpercii<=0; dpercii=0.5*10^-3; end 
        %d10ii=prctile(normrnd(d50ii,0.2*d50ii,1000,1),0.1)
        hydraulicData(ii,3)=dpercii;
        for jj=[find(isfinite(II(ii,:))==1)] % loop through all reaches downstream of reach ii

            reach_jj=find(raw_data_sort(:,ID_FromN)==jj); % position of reach j in the data
            if numel(reach_jj)>1
            reach_jj=reach_jj(1); %trough the resorting, some reaches are duplicate, use only the first value
            end

            
            if isempty(reach_jj)==0 % check if jj exists
              
                Dmat(ii,jj)=dii(ii);
                D90mat(ii,jj)=d90ii(ii);
                Thetamat(ii,jj)=2*d90ii(ii);
                Wmat(ii,jj)=raw_data_sort(reach_jj,ID_Wac);
                Lmat(ii,jj)=raw_data_sort(reach_jj,ID_Length);
            end % check if jj exists
        end  % loop through all reaches downstream of reach ii 

end % loop through all observations 

delete(h)