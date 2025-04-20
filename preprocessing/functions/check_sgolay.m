% double check the smoothing effect

clear
close all
clc

% pick a random participant/run
subj = 'p18_1';
run = 'sMGS1';
data_dir = fullfile('/mnt/Data1/eye_data/vswmda_eye/v3_short',run,subj,'result');
cd(data_dir)
file = dir('*step6.mat');
load(file.name)

% Lowpass filter window length
smoothInt = 20/1000; % in seconds

% Span of filter
span = ceil(smoothInt*edf.record.sample_rate);

% Construct the filter
%--------------------------------------------------------------------------
N = 2;                 % Order of polynomial fit
F = 2*ceil(span)-1;    % Window length
% [b,g] = sgolay(N,F);   % Calculate S-G coefficients

x = edf.samples.x_orig(edf.samples.trial == 1,set.eye);
y = edf.samples.y_orig(edf.samples.trial == 1,set.eye);

% X1 = conv(x,g(:,1), 'same');
% Y1 = conv(y,g(:,1), 'same');
X = smoothdata(x,'sgolay',F,'includenan','Degree',N);
Y = smoothdata(y,'sgolay',F,'includenan','Degree',N);

% 
% X2 = conv(x,g(:,2), 'same');
% Y2 = conv(y,g(:,2), 'same');


%%
f = figure(1);clf
plot(x,'k','LineWidth',2)
hold on
plot(X,'r','LineWidth',2)
hold on
% plot(X2,'b')
legend({'original X','SG filtered X'})
box off


