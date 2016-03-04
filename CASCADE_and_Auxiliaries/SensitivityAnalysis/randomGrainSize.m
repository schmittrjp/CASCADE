function [ dout ] = RandomGrainSize( din, pd )
% RANDOMGRAINSIZE derive a vector of randomly changed grain sizes: 
% The dsiturbed grainsize is defined as ds_dist=ds*Theta

%%%% Inputs: 
% din: Original vector of grain sizes for all sources 

%%%% Outputs: 
% dout: Randomly changed vector of grain sizes 

%%% Functioning: 
% 1st: Create a vector of disturbance. The vector of disturbance is
% generated from a Pearson random distribution. The parameters of this
% distribution are again random. This results in a wide range of
% differently shaped disturbance distributions

%% statistical parameters of disturbance (quite sensitive)
meanD=mean(din); 
skewness=randn(1);
stddev=std(din);
kurt=randn(1)+3; if kurt<0; kurt=0; end; 

% 2nd: Normalize the distribution between a min and a max value. E.g. Theta
% should be between 1 and 5. This upper and lower bound can again be
% random. 

scale_min=0.1 ; % smallest final disturbance 
scale_max=rand(1)*2; % largest final disturbanc 

% start 
n=length(din); % number of elements in the grain size vector 

% create original Pearson distribution. 
dist_out=nan(n,1);

while sum(isnan(dist_out))>0 % check if the distribution parameters result in a meaningful distribution 
moments = {meanD,stddev,skewness,kurt};
[dist_out,type] = pearsrnd(moments{:},n,1); % create disturbance vector 
end 

%  Rescale Pearson distribution between min and max value 

dist_min=min(dist_out);
dist_max=max(dist_out);

dist_out= (dist_out - dist_min).*(scale_max-scale_min)./(dist_max-dist_min) + scale_min
% adapted from: http://stackoverflow.com/questions/13137348/scaling-range-of-values-with-negative-numbers

% hist(dist_out,100)

dout=dist_out.*din; % change grain size
%  hist(dout,100);
end

