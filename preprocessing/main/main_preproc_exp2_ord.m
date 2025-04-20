% Preprocessing data for experiment 1
% Experiment 1 contains v3_short and v3_short_delay
% By: Linjing Jiang
% Update on: 04/05/2024

% f114, s12_1, sMGS4, almost no data

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
data_dir = fullfile('./eye_data','v3_model_ord'); %v3_model_quad

% We want to process each data folder separately
file_dirs = dir(fullfile(data_dir,'sMGS*','s*')); 

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

for ff = [1:21,23:length(file_dirs)] %ff_ind%
    clearvars -except data_dir out_dir file_dirs ff ff_ind 
    close all
    file_dir = fullfile(file_dirs(ff).folder,file_dirs(ff).name); % file directories
    out_dir = fullfile(file_dir,'result'); % output directories 

    fprintf('\n');
    fprintf('ff %03.0f subject %s run %s starts...',ff,file_dirs(ff).name,erase(file_dirs(ff).folder,[data_dir,'/']));

%% Step 1: File Conversion (Edfmex conversion toolbox)
% Set up any subfolder under top_dir indicating different subjects or
% recording sessions, e.g.,'1000/'
% Make sure that there is a single edf file under this folder
%file_dir = 'pf4_1\'; % Here we do not have a subfolder
%file_dir = '';

% Import the .edf file if it has never been imported before
if isempty(dir(fullfile(out_dir,'*.mat'))) %[file_dirs(ff).name '.mat']
    continue;
%    edf1 = edf2mat(file_dir,out_dir);
else
    edf_name = dir(fullfile(out_dir,'*.mat'));
    load(fullfile(out_dir,edf_name(1).name),'edf1')
end

% After you've done this, you will go into the specific data folder with
% the edf file. At the same time, in the Workspace you will see a "EDF2MAT"
% object which contains all the eye data imported from the edf file

% Check the manual of EDF2MAT if you are interested in what variables of
% the "EDF2MAT" object mean

%% Step 2 Set up parameters and load eye data from the EDF2MAT object

% Set up all the analysis parameters using the "setting" script
% A script called "setting" will automatically open. Please change the
% analysis parameters directly in the "setting" script. After you finalize
% the change, close the script, click the command window and press any key
% to proceed.
% open setting
% pause
set = setting_exp1(edf1);
% you can change parameters here if you want
%set.noise.blink_extend = 200;

% Get basic recording parameters from the edf file: sampling rate, pupil type,
% record type, eye recorded
% Note that here we created a new structure called 'edf' for the first time and
% load some recording parameters from 'edf1' (the EDF2MAT object) to 'edf'
edf = get_params(edf1);

% Set up screen parameters
% Note that there are 7 inputs to the 'get_screen_size' function, including
% 'edf': stored eye data structure (Please don't change this!!!)
% '1' (use customized parameters) or '0' (use default screen parameters).
% If you use customed parameters, please enter the following in sequence:
% distance from the eyes to the screen (in cm)
% width (in cm), height (in cm), x resolution (in pixel), y resolution
% (in pixel) of the monitor
edf = get_screen_size(edf,1,set.screen.d,set.screen.w,set.screen.h, ...
    set.screen.xres,set.screen.yres);

% Then, we extract important eye data from the EDF2MAT object and copy
% them to 'edf'
edf = load_sample(edf1,edf,set,file_dir,out_dir);

% Smooth the gaze and pupil data using SG filter
edf = sgolay_smoothing(edf,set); 

% Calculate velocity and acceleration
edf = cal_velacc(edf,set);

% Finally, we get calibration results from 'edf1' to 'edf' using the
% "get_calib" function
edf = get_calib(edf1,edf);

% Load parameters and detect task epochs

% Load parameter file
edf = load_param_exp1_exp2(edf);

% Detect different trials and task epochs
edf = detect_epoch_exp1_exp2(edf,set);

% All done! Now save the raw data to a .mat file
% fprintf('\n');
% fprintf('saving .mat file...');

filename = file_dirs(ff).name; edf.ID = filename;
%save(fullfile(out_dir,[filename,'_step2']));

%% Step 3 Artifact detection
% Next, let's detect artifacts in the data. This step would yield
% "edf.trackloss", under which there are multiple fields:
% 1. 'blink_ind': index for detected blinks
% 2. 'missing_ind': index for missing data (including blinks)
% 3. 'outside_ind': index for gaze position out of the screen boundary
% (either horizontal or vertical)
% 4. 'ext_ind': Extremely large sizes of pupil
% 5. 'pvel_ind': Extremely large velocity of pupil
% 6. 'all_ind': a combination of all the artifacts above
% 7. 'perc': percentage of artifacts overall
% You can set up the definition of most of those artifacts in the setting
% script
% Also note that all these indexes are based on the 'edf.samples'. For
% instance, an index of 300 indicates the 300th. row (sample) in any of the
% edf.samples array

% Detect artifacts
edf = detect_artifact(edf,set);

% Plot artifacts (will generate and save data figures automatically in the
% output folder)
plot_artifact(edf,out_dir,set);

% Output trackloss data (will generate a .csv file containing all the
% trackloss data)
tbl = table(edf.trackloss.perc,edf.blink.num,'VariableNames',{'Percentage of sample with artifacts','Number of blinks'});
writetable(tbl,fullfile(out_dir,['trackloss_id',edf.ID,'.csv']));
clear tbl

% save the data
clearvars edf1 edf_file
% fprintf('\n');
% fprintf('saving .mat file...');
%save([data_dir,file_dir,oclear all

% % Now, you should visually inspect blinks (optional if you are analyzing gaze locations)
% % If you are doing blink analysis, this step is a MUST
% % You should manually check all the blinks and see if those are truly
% % blinks, not other types of artifacts
% miniEye_ver0;
% % However, at this point, you can only inspect the blinks (its onset and
% % offset) without manually editing it. In the future I will add more
% % editing function to the GUI so that you can modify the blink events
% % directly using the toolbox

%% Step 4 Data Cleaning

close all

% Then, we need to remove artifacts from both the gaze and the pupil data
edf = remove_artifact(edf,set);

% Plot the time courses after artifact removal(automatically generate
% figures)
plot_timecourse_clean(edf,out_dir,set);

% % Alternatively, you can also inspect the data with our GUI
% miniEye_ver0;

% % IF YOU ARE ANALYZING PUPIL SIZE:
% % Do spatial interpolation of the pupil data
% edf = do_interpolation(edf,set);

% % Then do baseline correction if needed
% % Here we gave an exmaple of baseline correction for the pupil size
% edf = baseline_correction(edf,set);

% % Now you can inspect all types of the data again. Compare different plots and see
% % what their main differences are!
% miniEye_ver0;

% save the data
% fprintf('\n');
% fprintf('saving .mat file...');
%save([data_dir,file_dir,out_dir,filename,'_step4']);

%% Step 6 First round of event segmentation

% detect saccades first
% Velocity-based saccade detection of saccade (IVT, Reviewed in Salvucci & Goldberg, 
% 2000; Erkelens & Vogels, 1995; Sen & Megaw, 1984)
% Saccades were identified by identifying periods with velocity in excess 
% of 30°/s AND acceleration in excess of 8000°/s^2 AND duration for at least 8 ms 
edf = detect_saccades(edf,set);

% then, detect fixations
% Dispersion-based fixation detection (I-DT, Reviewed in Salvucci & Goldberg, 
% 2000; Widdle 1984)
% Fixations were defined as periods that do not belong to a saccade AND
% with dispersion less than 2 dva AND
% duration longer than 100 ms. 
edf = detect_fixations(edf,set);

%% It's time to double check the preprocessed results
% miniEye_2024

%% Step 7 Drift correction
% perform drift correction
edf = drift_correction(edf,set);

% after drift correction, has to recalculate x/y position in dva, as well 
% as velocity and acceleration
edf = load_sample_after_dc(edf,set);

% after drift correction we have to classify the fixations and saccades
% again
edf = detect_saccades_after_dc(edf,set);
edf = detect_fixations_after_dc(edf,set);

% combine saccades that are too close to each other
edf = combine_saccades(edf,set);

% calculate endpoint fixation
edf = cal_saccades_endpoint(edf,set);

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