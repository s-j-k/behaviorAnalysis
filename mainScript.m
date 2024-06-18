% ##############
% # preprocess #
% ##############

% mice numbers to analyze

% where you want to save the data
% pathsave = 'P:\su\DATA\behaviorData\opto\cohort2\results';
% datapath = 'P:\su\DATA\behaviorData\opto\cohort2\';
pathsave = 'P:\su\DATA\behaviorData\opto\cohort1\results';
datapath = 'P:\su\DATA\behaviorData\opto\cohort1\';
% fill mouse metadata to excel spreadsheet
% settings={};

% data=preprocess();
% creates large data struct with fields


% #################
% # main analysis #
% #################

% load settings file from an excel spreadsheet
% this spreadsheet will be populated with the finished results
% results will also be saved to pathsave as .mat
% all info about mice (genotype, raw data path, etc.) will be loaded from the xlxs
% settings={};
% optoAnalysis(number_mice,pathsave,settings)

% mice numbers to analyze
number_mice = [37,39,41,42, 43, 44, 45, 46];
% number_mice=[49,50,51,52,53,54,55,56];
nFiles=23; %days of behavior
% nFiles=12;
plot_indiv_data = true;
plot_lick_psth = true;
savefig = true;
optoAnalysis(number_mice,pathsave,datapath,nFiles,plot_indiv_data,plot_lick_psth,savefig)

%%

% add a line to take out duplicate sessions from same day
% add a line to concatenate sessions together from same day 
% #################
% # lick analysis #
% #################

% takes data structure from preprocessing and gets lick raster plots
% takes inputs data (data struct), hidprime (1 or 0, to use only high d
% prime during probe, and pathsave for the figures
lickRaster(data,hidprime)

%%
%figures to keep
nums=1;
figs2keep = [nums];

all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));
% Uncomment the following to 
% include ALL windows, including those with hidden handles (e.g. GUIs)
% all_figs = findall(0, 'type', 'figure');