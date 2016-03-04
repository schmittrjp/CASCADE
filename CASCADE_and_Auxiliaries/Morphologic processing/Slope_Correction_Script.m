

%% slope correction for reaches with void measurements
% There are three methods

% Method 1 ('SSCorMeth 1'), adapats the slope of the void reach, and the
% slope of the downstream reaches, in order to derive a new slope for both.
%->SLOPE IS A Function of measured morphologic conditions, but break points
%are smoothened. 

% Method 2 ('SSCorMeth 2'), uses a fixed increase in the upstream
% elevation, to locally create a very low slope. 

% Method 3: ('SSCorMeth 3'), uses a correlation with the drainage area in order ot derive local slopes. 
       switch SSCorMeth
  
        case  'SSCorMeth 3'
                       
           ss_val=find(raw_data_sort(:,ID_Slp)>1e-5);% find valid slope indices
           x=raw_data_sort(ss_val,ID_Ad); 
           y=raw_data_sort(ss_val,ID_Slp);

           % Set up fittype and options.
           ft = fittype( 'power1' ); opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); opts.Display = 'Off';
           opts.StartPoint = [0.0528236677872659 -0.169343731781852];

           % Fit model to data.
           [fitresult_slope, gof] = fit( x, y, ft, opts );

       end

%% begin the correction 
SS_cor_results=zeros(length(raw_data_sort),3);
ss=1;
for sss=1:size(raw_data_sort,1)  %% loop through all reaches 

SS_cor_results(ss,1:2)=[sss raw_data_sort(sss,ID_Slp)]; 

if raw_data_sort(sss,ID_Slp)<=1e-5 % find reaches with a invalid slope
ss=ss+1;
    switch SSCorMeth
  
    case  'SSCorMeth 1'
         waitbar(ss/sum(raw_data_sort(:,ID_Slp)<1e-5)); 
 
        
           % Find the path from the current node sss to the outlet
           [~,s_dsnodes]=graphshortestpath(D,sss,48865);
        
           s_us_elev =raw_data_sort(s_dsnodes,ID_ElDs) ; % get the upstream elevation for that downstream nodes
           delta_h=(s_us_elev-raw_data_sort(sss,ID_ElDs)); % elevation difference between upstream and downstrem nodes
           c_correction_nodes=s_dsnodes(delta_h<-1); % find all downstream nodes that have a lower elevation; 

           if isempty(c_correction_nodes)==0
           c_correction_nodes=c_correction_nodes(1); % find the first node on the downstream path that has a lower slope 
           [~,s_dsnodes]=graphshortestpath(D,sss,c_correction_nodes);
           
           %calculate new slope
           sss_new_slp=(raw_data_sort(sss,ID_ElUs)-raw_data_sort(c_correction_nodes(end),ID_ElDs))/... % new sloep is the elevation distance between the two reaches  
          graphshortestpath(D,sss,c_correction_nodes); % divided by the combined length

            if sss_new_slp<0
                kli=3
            end
      
           %%%% method 1 
           % this smoothes out all breaks: 
           raw_data_sort([sss s_dsnodes(1:end-1)],ID_Slp)=sss_new_slp; % replace the slope of all downstream nodes required for correction

           end
        
    case  'SSCorMeth 2'

           %%%% method 2 
           % this is not based in the data, 
           raw_data_sort(sss,ID_Slp)=0.25/raw_data_sort(sss,ID_Length);     
            
    case  'SSCorMeth 3'

           %%%% method 2 
           % this is not based on a extrapolatiioon from Q1.5
           raw_data_sort(sss,ID_Slp)=fitresult_slope(raw_data_sort(sss,ID_Q15));    
    

    end 
SS_cor_results(ss,3)=[raw_data_sort(sss,ID_Slp)];
  
end
end

SS_cor_results(ss:end,:)=[];
