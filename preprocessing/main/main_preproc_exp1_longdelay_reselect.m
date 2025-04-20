% Preprocessing data for experiment 1
% Experiment 1 contains v3_short and v3_short_delay
% By: Linjing Jiang
% Update on: 04/05/2024

%% Step 0: Preparation
% The first step is to import the edf file into MATLAB, using the
% Edf2Mat toolbox
% toolbox

% Clean the workspace and everything
clear all
close all
clc

% v3_short:
% Set up the data directory, where you store all your data,
% e.g., 'D:/linjing_eyetracking/test'
data_dir = fullfile('./eye_data/','v3_short'); %v3_model_quad

% We want to process each data folder separately
file_dirs = dir(fullfile(data_dir,'sMGS*','p*')); 


addpath(genpath(data_dir))

% Set up the script directory, where you store the scripts
script_dir = fullfile('./scripts','Cloned-Repo');
addpath(genpath(script_dir))

% f215, s43, sMGS6, trial 10 - saccade has lower acceleration than
% threshold
% f298, s54, sMGS8
% f7, p18, sMGS1
% f13, p24, sMGS1
% f136, p38, sMGS4
% f92, p30, sMGS3

% ff 92 - file corrupted, need to rerun
% ff 100 - file corrupted, need to rerun
% ff 251 - file corrupted
for ff = 1:length(file_dirs)%[1:91,93:99,101:length(file_dirs)] %ff_ind%
    clearvars -except data_dir out_dir file_dirs ff ff_ind 
    close all
    file_dir = fullfile(file_dirs(ff).folder,file_dirs(ff).name); % file directories
    out_dir = fullfile(file_dir,'result'); % output directories 

    fprintf('\n');
    fprintf('ff %03.0f subject %s run %s starts...',ff,file_dirs(ff).name,erase(file_dirs(ff).folder,[data_dir,'/']));


% Import the .edf file if it has never been imported before
if isempty(dir(fullfile(out_dir,'*step6.mat'))) %[file_dirs(ff).name '.mat']
    continue;
%    edf1 = edf2mat(file_dir,out_dir);
else
    edf_name = dir(fullfile(out_dir,'*step6.mat'));
    load(fullfile(out_dir,edf_name(1).name),'edf','set')
end

filename = edf.ID;
%% Step 8 Select events of interest
% Select events of interests
edf = select_saccades_exp1_exp2(edf,set); % select saccades

%% STep 9 Further cleaning
% Detect trials breaking fixation
edf = detect_breakfix(edf);

% Calculate angular and order error
edf = cal_order_err(edf);

% save the data
% fprintf('\n');
% fprintf('saving .mat file...');
save(fullfile(out_dir,[filename,'_step6']));

% % Timing
% timing = [nan(1,13);edf.timing.epoch_dur_diff;nan(1,13)];
% timing = [timing edf.timing.tr_dur_diff];
% writematrix(timing,'timing.csv')

end