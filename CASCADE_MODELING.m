%%%% CASCADE Master script. This script contains all necessary sub-scripts
%%%% and functions to either 
% a) Run the model from the raw data
% b) Load data from my previous experiments and look only at specific parts
% of the modeling process
% c) Make runs of the CASCADE model and look at its results. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define output folder for the preprocessed data
PreprocessedOutFolder='CASCADE_PreProcessing_Output'; 

% This was added as a test for GitHub
% hello hello
% this is new stuff

addpath(genpath(pwd))

% Run the master pre-processing from raw data, or load preprocessed data from a file 
    prompt = 'Do you want to rerun preprocessing? [Y/N]: ';
    str = input(prompt,'s');

    if strcmpi('y',str)
        
  