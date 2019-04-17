%% call get_MS

clear all
close all
clc

% load data
load('/data/home/ebalestr/PhD/my_toolboxes/sandbox_data2.mat')

% define configuration
cfg.screensize = [1920 1080]; 
cfg.pix_per_va = 51;
cfg.kernellength = 3;
cfg.toi = [-2 0];

MS_data = get_MS(s_data, cfg);

lgcl_arb_answ_RIGHT = MS_data.trialinfo(:,3)==55 & MS_data.trialinfo(:,2)==88 ...
    & MS_data.trialinfo(:,1) == 19;
lgcl_arb_answ_LEFT = MS_data.trialinfo(:,3)==55 & MS_data.trialinfo(:,2)==88 ...
    & MS_data.trialinfo(:,1) == 91;

cfg_s.trial = lgcl_arb_answ_LEFT;
reduced_data_LEFT = trialselect(MS_data, cfg_s);

cfg_s.trial = lgcl_arb_answ_RIGHT;
reduced_data_RIGHT = trialselect(MS_data, cfg_s);

angles_left = reduced_data_LEFT.angle(reduced_data_LEFT.lgcl_MS_onset);
angles_right = reduced_data_RIGHT.angle(reduced_data_RIGHT.lgcl_MS_onset);

angles_left = reduced_data_LEFT.angle(reduced_data_LEFT.lgcl_mask_MS);
angles_right = reduced_data_RIGHT.angle(reduced_data_RIGHT.lgcl_mask_MS);

figure; 
polarhistogram(angles_right); hold on
polarhistogram(angles_left)

avg_left = rad2deg(circmean(angles_left));
avg_right = rad2deg(circmean(angles_right));

