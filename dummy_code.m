clear
clc
% close all


addpath('/home/elio/toolboxes/eeglab2021.1')
eeglab nogui

% get filelist (in case of loop)
filelist = strsplit(ls('/home/elio/PhD/alphaReturn/fulldata/*.set'));
filelist = sort(filelist); 
filelist(1) = []; % linux specific(?)

%% cfg definition & MS computation

% MS
cfg=[];
cfg.kernellength = 3; % minimal MS duration (in data samples)
cfg.noisesuppress = true; % default true, check help
cfg.lambda = 5; % default 5, check help
cfg.labelMS = {[-1, -.1], [.1, 1]}; % cosine value for left vs right definition, check help
cfg.toi = [-200, 1500]; % time of interest

% MS rate
cfg_MSrate = [];
cfg_MSrate.fsample = 256;
cfg_MSrate.movingwin = .1;

%% apply MS function


% load EEG dat
nsubjs = length(filelist); isubj = 16;
fname = filelist{isubj};
EEG = pop_loadset(fname);

% compute MS
[MSstruct, s_data] = get_MS(EEG, cfg);

% define experimental condtions (based on trialinfo)
% obviously experiment-specific 
table_tinfo = struct2table(EEG.trialinfo);
trial_SOA = table_tinfo.soa;
unSOA = unique(trial_SOA); %unSOA = unSOA(2:end);
longSOAs_mask = (trial_SOA == 0);
left_cue_mask = table_tinfo.cue_position == 1;
right_cue_mask = table_tinfo.cue_position == 2;

% definition of logical masks to be used for assigning trials to one
% condition or the other
mask_long_left = longSOAs_mask & left_cue_mask;
mask_long_right = longSOAs_mask & right_cue_mask;

% compute MS rate for both conditions
MSrate_L = get_MSrate_EXPcondition(s_data, mask_long_left, cfg_MSrate);
MSrate_R = get_MSrate_EXPcondition(s_data, mask_long_right, cfg_MSrate);

plot_summary(s_data)



