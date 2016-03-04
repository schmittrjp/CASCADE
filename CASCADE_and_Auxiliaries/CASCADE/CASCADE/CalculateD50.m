function [d50] = calculateD50(dmat,flowmat)
%CALCULATED50 Calculates the median grain size of a reach. 

% Inputs:
% dmat: matrix of grain sizes that can potentially be found in a reach. 
% flowmat: matrix of flows (e.g. inputs, outputs, or sediment balance) by
% which to weigh the observations in the d50 calculation. 

% Outputs: 
% d50: vector with d50 for all reaches. 

weightMat=flowmat./nansum(flowmat); 

weightedMedian(dmat,weightMat)




end

