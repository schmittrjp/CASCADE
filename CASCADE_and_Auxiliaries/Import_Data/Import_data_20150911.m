%%%%% Import data 20150817 

%% First import the raw data matrix 

% Import the data
[~, ~, raw] = xlsread('\CASCADE_and_Auxiliaries\Import_Data\20150911_Final_River_Net.xlsx','Data','A2:T49806');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
raw_data = reshape([raw{:}],size(raw));

% Clear temporary variables
clearvars raw R;

%% Second, import Name and columns (ID) of each variables 


% Import the data
[~, col_names, raw] = xlsread('E:\Connectivity classification Da\Matlab\Import_Data\20150911_Final_River_Net.xlsx','Variable Names');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));
data(1,:)=[];
% Allocate imported array to column variable names
global ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY ID_SubWS % Clear temporary variables


ID_arcid = data(:,1);
ID_FromN = data(:,2);
ID_ToN = data(:,3);
ID_ElUs = data(:,4);
ID_ElUsRaw = data(:,5);
ID_ElDs = data(:,6);
ID_ElDsRaw = data(:,7);
ID_Slp = data(:,8);
ID_SlpRaw = data(:,9);
ID_ElDiff = data(:,10);
ID_Length = data(:,11);
ID_StrO = data(:,12);
ID_MicroWSAre = data(:,13);
ID_FldPlnWdth = data(:,14);
ID_Ad = data(:,15);
ID_FX = data(:,16);
ID_FY = data(:,17);
ID_TX = data(:,18);
ID_TY = data(:,19);
ID_SubWS = data(:,20);
% clearvars data raw R;