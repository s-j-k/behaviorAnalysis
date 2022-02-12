%% preprocess

% mice numbers to analyze

% where you want to save the data

% fill mouse metadata to excel spreadsheet
settings={};

%% run initial analysis

% mice numbers to analyze
number_mice = [37,39,41,42, 43, 44, 45, 46];

% where you want to save the data
pathsave = 'T:\su\DATA\behaviorData\opto\cohort1\results';

% load settings file from an excel spreadsheet
% this spreadsheet will be populated with the finished results
% results will also be saved to pathsave as .mat
% all info about mice (genotype, raw data path, etc.) will be loaded from the xlxs
settings={};
optoAnalysis(number_mice,pathsave,settings)

%%
