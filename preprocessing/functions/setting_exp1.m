function set = setting_exp1(edf)
% settings for eye-tracking analysis
% please change the value directly below based on your analysis preferences

%% The sequence of the messages within a trial
set.msg = {'TRIALID','FixSrt','Tar1Srt','Gap1Srt','Tar2Srt','Gap2Srt',...
    'Tar3Srt','Gap3Srt','Tar4Srt','DelaySrt','CueSrt','SacSrt','ItiSrt','TRIAL_RESULT 0'};
% In the following analysis, each of these messages will be marked as from
% 1 to n

%% Which eye to analyze
% 1. left
% 2. right
% 3. binocular

% Check Left
le = isequal(edf.Samples.gx(:,1),-32768*ones(size(edf.Samples.gx(:,1))));
% Check Right
re = isequal(edf.Samples.gx(:,2),-32768*ones(size(edf.Samples.gx(:,2))));

if le & ~re
set.eye = 2; 
elseif ~le & re
    set.eye = 1;
elseif ~le & ~re
    set.eye = 2; % default, right eye
elseif le & re
    error('No valid data!');
end

%% Screen info
% distance from the eyes to the screen (in cm)
set.screen.d = 67.7;

% width (in cm)
set.screen.w = 38;

% height (in cm)
set.screen.h = 30.5;

% x resolution (in pixel)
set.screen.xres = 1280;

% y resolution (in pixel) of the monitor
set.screen.yres = 1024;

% % fMRI
% % distance from the eyes to the screen (in cm)
% set.screen.d = 143.5;
% 
% % width (in cm)
% set.screen.w = 43;
% 
% % height (in cm)
% set.screen.h = 35;
% 
% % x resolution (in pixel)
% set.screen.xres = 1024;
% 
% % y resolution (in pixel) of the monitor
% set.screen.yres = 768;

%% Saccade and fixation detection methods
% % Nystrom's detection algorithm
% % minimum saccade duration
% set.sac.min_dur = 0.15; % 150 ms
% 
% % initial peak velocity threshold
% set.sac.init_velpeak_threshold = 100; % Initial value of the peak detection threshold. 
% 
% % minimum fixation duration
% set.fix.min_dur = 0.030; % in seconds
% 
% % blink velocity threshold
% set.blink.vel_threshold = 1000; % if vel > 1000 degrees/s, it is noise or blinks
% 
% % blink acceleration threshold
% set.blink.acc_threshold = 100000; % if acc > 100000 degrees/s^2, it is noise or blinks

% Velocity-based saccade detection (IVT, Reviewed in Salvucci & Goldberg, 
% 2000; Erkelens & Vogels, 1995; Sen & Megaw, 1984)
% Saccades were identified by identifying periods with velocity in excess 
% of 30°/s and acceleration in excess of 8000°/s^2 for at least 8 ms and 
% resulting in at least a 0.25° amplitude gaze shift.
set.sac.vel_threshold = 30; %threshold velocity for saccades (deg/s)
set.sac.acc_threshold = 6000; %threshold acceleration for saccades (deg/s^2)
set.sac.amp_threshold = 0.25; %threshold amplitude for saccades (deg)
set.sac.dur_threshold = 8; % duration threshold for saccades (ms)
set.sac.comb_threshold = 50; % combine saccades that are too close to each other (current offset and next onset is less than 25 ms apart)
% How to combine:
% the first saccade startpoint would be the new saccade startpoint
% the last saccade endpoint would be the new saccade's endpoint

% Saccade selection criteria
set.sac.srt_err = 2; % we want to select saccades with its start position less than 2 dva from the center cross
set.sac.end_err = 3.5; % we want to select saccades with its end position less than 3.5 dva from the target or nontargets
set.sac.aend_err = 2; % the closest saccade's end position whould be more than 2 dva from the center cross
set.sac.msg = [11 12]; % we want to select saccades that start after message 11 or 12 (cue or saccade)
set.sac.rt = [50 1200]; % primary saccade latency should be longer than 50 ms and shorter than 1200 ms
% did not use the following criteria
% set.sac.fixdur = 50; % we want to select saccades followed by a fixation of at least 50 ms

% Dispersion-based fixation detection (I-DT, Reviewed in Salvucci & Goldberg, 
% 2000; Widdle 1984)
% fixations are defined as samples with a velocity lower than or equal to 
% 10 dva/s, duration longer than or equal to 30 ms, and both horizontal and 
% vertical gaze dispersion smaller than or equal to 1.5 dva
set.fix.disp_threshold = 1.5; % dispersion threshold for fixations (dva)
set.fix.dur_threshold = 30; % duration threshold for fixations (ms)
% Have a lower velocity threshold to define fixations
set.fix.vel_threshold = 10;

%% Blink detection methods

% set.blink.methods
% methods of performing blink detection

% 1. default blink detection from eyelink
% Definition: a period of saccade detector activity with the pupil data 
% missing for three or more samples in a sequence

% 2. default blink detection from eyelink with extension of time window
% set.blink.extend: the time window extended before and after the blink (in
% ms)

% 3. modified velocity algorithm by Nystrom 
% Nyström, M., & Holmqvist, K. (2010). An adaptive algorithm for fixation, 
% saccade, and glissade detection in eyetracking data. Behavior research methods, 42(1), 188-204.

% 4. noise-based algorithm by Hershman et al. (2018)

set.noise.blink_method = 3; 
set.noise.blink_extend = 100; % 200 ms before and after detected blink. Will only affect method 2, 3 and 4. (100ms is good for method 4).
set.noise.blink_pvel = 8000; % pupil velocity threshold (> 5000 would be identified as a blink). Will only affect method 3.

set.noise.gaze_vel = 1000; % if gaze velocity > 1000 degrees/s, it is noise or blinks
set.noise.gaze_acc = 100000; % if gaze acceleration > 100000 degrees/s^2, it is noise or blinks
% see Duchowski, 2003, for an overview on velocity and acceleration thresholds

set.noise.pupil_sz = 100000; % how many median absolute deviations of pupil size - not used
set.noise.pupil_vel = 100000; % how many median absolute deviations for pupil velocity - not used

%set.extreme_pupil.lamda = 3; % how many median absolute deviations for pupil size
%set.extreme_pupil_vel.lamda = 100; % how many median absolute deviations for pupil change speed
%set.extreme_pupil.threshold = 1000;             
%ETparams.blinkAccThreshold = 100000;               

%% Drift correction threshold
%set.noise.dc_threshold = 2; % if fixation position is more than 2 dva from the center cross, we no longer perform drift correction
set.noise.baseline_ratio = 1/2; % baseline is half of the fixation duration right before the stimulus onset 
% value could range from 0 to 1 (could not be 0).

%% Data cleaning options
set.clean.artifact = 0; % choose which type of artifact to remove from the data
% 0: all artifacts
% 1: blinks
% 2: missing data
% 3. gaze data outside of the screen boundary
% 4. extreme velocity and acceleration of the gaze
% 5. extreme size of the pupil
% 6. extreme velocity of the pupil

set.clean.dist = 50; % if two artifacts are 50 ms within each other, then merge these two 
%% Interpolation (for pupil data)
set.interp.method = 1; 
% 1: linear interpolation
% 2: spine interpolation

set.interp.range = 50; % Interpolation based on how many ms samples prior and after

%% Baseline correction (for pupil data)
set.bcorr.msg = 2; % choose the period between the 2nd and the 3rd message as the baseline period
% In this scenario, it is the fixation period

% choose the method of baseline correction
set.bcorr.method = 1; % 1, substractive; 2, divisive

% choose how to calculate the baseline pupil value, by taking either mean
% or median of the baseline period throughout the task
set.bcorr.base_method = 2; % 1, mean; 2, median of the baseline values

%% Which trial to plot
set.plot.trial_srt = 1; % start trial of plotting
set.plot.trial_end = 5; % end trial of plotting

end

