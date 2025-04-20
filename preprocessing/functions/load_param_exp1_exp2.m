function edf = load_param_exp1_exp2(edf)
% Load parameter files
% Make sure you modify this function based on different experiments
% Go to the specific data folder
data_f = dir(fullfile(edf.dir.file_dir,'*TRIAL*.dat'));
if isempty(data_f)
    data_f = dir(fullfile(edf.dir.file_dir,'*Trial*.dat'));
end

params = [];
for ii = 1:length(data_f)
    if verLessThan('matlab','9.8')
            % for r2019b or earlier version
params = [params;readtable(fullfile(data_f(ii).folder,data_f(ii).name),'PreserveVariableNames',1)];
    else
    % for newer versions
params = [params;readtable(fullfile(data_f(ii).folder,data_f(ii).name),'VariableNamingRule','preserve')];
    end
end

% 26:33: x and y location
% 4-7: quadrant order
% 46:57: timing

% 34: probe order
% 35: probe quadrant
% 36: probe image category
% 42,43: probe x and y location

x = params(:,26:29);
y = params(:,30:33);
ord = params(:,4:7);
tm = params(:,[46:57]);

tarx = params(:,42);
tary = params(:,43);
tar_ord = params(:,34);
tar_quad = params(:,35);
tar_img = params(:,36);
ecc = params(:,65);

% Store all the info
edf.param.datasource = (params);
edf.param.stimx = table2array(x);
edf.param.stimy = table2array(y); edf.param.stimy = edf.screen.yres - edf.param.stimy; % flip y
edf.param.tarx = table2array(tarx);
edf.param.tary = table2array(tary); edf.param.tary = edf.screen.yres - edf.param.tary; % flip y
[edf.param.stimx_deg,edf.param.stimy_deg] = pix2ang(edf.param.stimx,edf.param.stimy,edf);
[edf.param.tarx_deg,edf.param.tary_deg] = pix2ang(edf.param.tarx,edf.param.tary,edf);
edf.param.time = table2array(tm)*1000; % change to miliseconds
edf.param.quad_order = table2array(ord);
edf.param.probe_order = table2array(tar_ord);
edf.param.probe_quad = table2array(tar_quad);
edf.param.probe_img = table2array(tar_img);
% eccentricity
            edf.param.ecc = table2array(ecc);

end