% Figure 1 for serial order behavioral paper
% By: Linjing Jiang
% Date: 06/10/2024

%% Fig. 1C plot target locations, spatial and temporal patterns
clear 
close all
clc

% target locations - behavioral experiments
ang = [15:20:75];%11.25:22.5:360;
ang = [ang ang+90 ang+180 ang+270];
ecc = [3 5.5];
x = [ecc(1)*cosd(ang) ecc(2)*cosd(ang)];
y = [ecc(1)*sind(ang) ecc(2)*sind(ang)];
figure(1);clf
scatter(x,y,50,'filled','k')
hold on
line([0 0],[-6 6],'Color','k')
line([-6 6],[0 0],'Color','k')
axis square
box off
axis off
print(figure(1),'Fig1_loc.png','-dpng','-r300')

%% Fig. 1C plot spatial and temporal patterns

%d = 138.9; % cm
%h = 35;
%w = 43;
%x_res = 1024;
%y_res = 768;

d = 67.7; % cm
h = 30.5;
w = 38;
x_res = 1280;
y_res = 1024;


% construct different combinations
data = nchoosek(ang,4);

% get rid of angles in the same quadrant
quad = zeros(size(data,1),1);
for ii = 1:size(data,1)
    v = quadrant(data(ii,:));
    if length(v) == length(unique(v)) % if in different quadrants
        quad(ii) = 0;
    else
        quad(ii) = 1;
    end
end
data_quad = data(quad==0,:);

% include fully asymmetric + only one axial symmetric + only one central symmetric case
s = zeros(size(data_quad,1),1);
for ii = 1:size(data_quad,1)
    s(ii) = symmetric5(data_quad(ii,:));
end
data_quad_sym = data_quad(s==0,:);

% If you want to check symmetric type, uncommend the following section
% check symmetric type
% include fully asymmetric + only one axial symmetric + only one central symmetric case
st = zeros(size(data_quad,1),1);
for ii = 1:size(data_quad,1)
    st(ii) = symmetric6(data_quad(ii,:));
end
% plot different types of location combinations
% st == 1: axial symmetric
% st == 2: central symmetry
% st == 3: fully assymetric
ind = find(st==1); ind = ind(1);
x = ecc(2)*cosd(data_quad(ind,:));
y = ecc(2)*sind(data_quad(ind,:));
X = ecc(2)*cosd(ang);
Y = ecc(2)*sind(ang);
figure(1);clf
scatter(X,Y,600,'k','LineWidth',1)
hold on
scatter(x,y,600,'filled','k','LineWidth',3)
xlim([-7 7])
ylim([-7 7])
line([0 0],[-7 7],'Color','k','LineStyle','--','LineWidth',3)
line([-7 7],[0 0],'Color','k','LineStyle','--','LineWidth',3);
set(gca,'PlotBoxAspectRatio',[w h 1])
axis off
axis square
print(figure(1),'Fig1_AS.png','-dpng','-r300')

ind = find(st==2); ind = ind(17);
x = ecc(2)*cosd(data_quad(ind,:));
y = ecc(2)*sind(data_quad(ind,:));
X = ecc(2)*cosd(ang);
Y = ecc(2)*sind(ang);
figure(1);clf
scatter(X,Y,600,'k','LineWidth',1)
hold on
scatter(x,y,600,'filled','k','LineWidth',3)
xlim([-7 7])
ylim([-7 7])
line([0 0],[-7 7],'Color','k','LineStyle','--','LineWidth',3)
line([-7 7],[0 0],'Color','k','LineStyle','--','LineWidth',3);
set(gca,'PlotBoxAspectRatio',[w h 1])
axis off
axis square
print(figure(1),'Fig1_CS.png','-dpng','-r300')

ind = find(st==3); ind = ind(20);
x = ecc(2)*cosd(data_quad(ind,:));
y = ecc(2)*sind(data_quad(ind,:));
X = ecc(2)*cosd(ang);
Y = ecc(2)*sind(ang);
figure(1);clf
scatter(X,Y,600,'k','LineWidth',1)
hold on
scatter(x,y,600,'filled','k','LineWidth',3)
xlim([-7 7])
ylim([-7 7])
line([0 0],[-7 7],'Color','k','LineStyle','--','LineWidth',3)
line([-7 7],[0 0],'Color','k','LineStyle','--','LineWidth',3);
set(gca,'PlotBoxAspectRatio',[w h 1])
axis off
axis square
print(figure(1),'Fig1_FS.png','-dpng','-r300')

%% Fig. 1D plot eye movement trajectory
% run plot_single_trial.m
clear
close all
clc

% pick a random participant/run
subj = 'p13_1'; % p18_1
run = 'sMGS6';
data_dir = fullfile('./v3_short',run,subj,'result');
file = dir(fullfile(data_dir,'*step6.mat'));
load(fullfile(file.folder,file.name))

trial_property = 8; % trial 13
msg_property = 'On';
sac_property = 'Off';
fix_property = 'Off';

time = (edf.samples.time - edf.samples.time(1))/1000; % time

p = edf.samples.pupil_size(:,set.eye); % pupil
X = edf.samples.x_deg(:,set.eye); % Gaze X
Y = edf.samples.y_deg(:,set.eye); % Gaze Y
% X = edf.samples.x(:,set.eye); % Gaze X
% Y = edf.samples.y(:,set.eye); % Gaze Y

p_clean = edf.samples.pupil_size_clean(:,set.eye); % pupil
% X_clean = edf.samples.x_deg_clean(:,set.eye); % Gaze X
% Y_clean = edf.samples.y_deg_clean(:,set.eye); % Gaze Y

X_clean_dc = edf.samples.x_deg_clean_dc(:,set.eye); % Gaze X
Y_clean_dc = edf.samples.y_deg_clean_dc(:,set.eye); % Gaze Y
% X_clean_dc = edf.samples.x_clean_dc(:,set.eye); % Gaze X
% Y_clean_dc = edf.samples.y_clean_dc(:,set.eye); % Gaze Y

tarX = edf.param.tarx_deg(trial_property);
tarY = edf.param.tary_deg(trial_property);
stimX = edf.param.stimx_deg(trial_property,:);
stimY = edf.param.stimy_deg(trial_property,:);
taridx = edf.param.probe_order(trial_property);
% tarX = edf.param.tarx(trial_property);
% tarY = edf.param.tary(trial_property);

% index of a trial
ind_trial = find(edf.samples.trial == trial_property);
% start time of a trial
srt_trial = time(ind_trial(1));
end_trial = time(ind_trial(end));

%srt_property = 0; % we want to set the trial starting time to zero
%end_property = time(ind_trial(end))-srt_trial; %time(ind_trial(1));

% start time of response period
srt_resp = (edf.events.msg.time(trial_property,11) - edf.samples.time(1))/1000

% plot specific time window within a trial
t1 = srt_trial;%srt_property + srt_trial; % when do this trials actually start?
t2 = end_trial;%end_property + srt_trial; % when do this trial actually end?
t3 = srt_resp;% when the response period starts?

% find index of those times
ind = find(time >= t1 & time <= t2);

% find index of all the messages
n_msg = max(edf.samples.msg);
%             for ii = 1:n_msg
%                 indm(ii,:) = find(edf.samples.msg == ii)';
%             end
indm = edf.events.msg.ind_srt;

% plot
f = figure('Renderer', 'painters', 'Position', [10 10 1500 400])
% gaze x & y
% plot target locations

%line([0 t2-t1],[stimX(taridx) stimX(taridx)],'LineStyle','--','Color',[0, 0.4470, 0.7410,1])
%hold on
%line([0 t2-t1],[stimY(taridx) stimY(taridx)],'LineStyle','--','Color',[0.8500, 0.3250, 0.0980,1])
%hold on
p1 = line([t3-t1 t2-t1],[tarX tarX],'LineStyle','--','Color', ...
    [0, 0.4470, 0.7410,1],'LineWidth',2)
hold on
p2 = line([t3-t1 t2-t1],[tarY tarY],'LineStyle','--','Color', ...
    [0.8500, 0.3250, 0.0980,1],'LineWidth',2)
hold on

% minY and maxY
minY = min([X(ind);Y(ind);tarX;tarY;X_clean_dc(ind);Y_clean_dc(ind)])-0.5;
maxY = max([X(ind);Y(ind);tarX;tarY;X_clean_dc(ind);Y_clean_dc(ind)])+0.5;
% gaze locations
p3 = plot(time(ind)-t1,X(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,0.2]);
hold on
p4 = plot(time(ind)-t1,Y(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,0.2]);
% plot(time(ind)-t1,X_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
% plot(time(ind)-t1,Y_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
p5 = plot(time(ind)-t1,X_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,1]);
p6 = plot(time(ind)-t1,Y_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,1]);
xlim([0 t2-t1])
%ylim([-4 4])
%ylim([minY maxY])
xlabel('Time (s)')
ylabel('Eye Position (dva)')
legend('off');
%title(['Gaze time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
hold on

% plot messages
minY = -4;
maxY = 5;
if strcmp(msg_property,'On')
    for ii = trial_property%1:size(indm,1)
        for jj = [3,10,11]%1:size(indm,2)
            h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[minY,maxY], ...
                'LineStyle','--','Color','k','LineWidth',2);% cmap(ii,:)
        %    text(time(indm(ii,jj))-t1,maxY*0.9,0,num2str(edf.events.msg.txt(ii,jj)));
        end
    end
end
hold on

% plot saccades

% plot only primary saccade
ind_psac = edf.cal.primary_sac_ind(trial_property);
if ~isnan(ind_psac)
    psac_on = edf.events.sac_dc.ind_srt(ind_psac);
    psac_off = edf.events.sac_dc.ind_end(ind_psac);
 %   plot(time(psac_on:psac_off)-t1,X_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
 %   plot(time(psac_on:psac_off)-t1,Y_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
% endpoint fixations
    pfix_on = edf.events.fix_dc.ind_srt(edf.events.sac_dc.endfix(ind_psac));
    pfix_off = edf.events.fix_dc.ind_end(edf.events.sac_dc.endfix(ind_psac));
p7 = plot(time(pfix_on:pfix_off)-t1,X_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
p8 = plot(time(pfix_on:pfix_off)-t1,Y_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
end

% plot most accurate saccade
ind_psac = edf.cal.acc_sac_tar_ind(trial_property);
if ~isnan(ind_psac)
    psac_on = edf.events.sac_dc.ind_srt(ind_psac);
    psac_off = edf.events.sac_dc.ind_end(ind_psac);
%    plot(time(psac_on:psac_off)-t1,X_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
%    plot(time(psac_on:psac_off)-t1,Y_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);

% endpoint fixations
    pfix_on = edf.events.fix_dc.ind_srt(edf.events.sac_dc.endfix(ind_psac));
    pfix_off = edf.events.fix_dc.ind_end(edf.events.sac_dc.endfix(ind_psac));
    plot(time(pfix_on:pfix_off)-t1,X_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
    plot(time(pfix_on:pfix_off)-t1,Y_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);

end


% plot fixations used to perform drift correction
idx_fix = find(edf.samples.baseline_dc == 1 & ... % is a baseline
    edf.samples.trial == trial_property); 

dif_fix = diff(idx_fix);
idx_fix_sep = find(dif_fix > 1); idx_fix_sep = [0;idx_fix_sep;length(idx_fix)];
for xx = 1:(length(idx_fix_sep)-1)
    idx_temp = (idx_fix_sep(xx)+1):(idx_fix_sep(xx+1));
    plot(time(idx_fix(idx_temp))-t1,X_clean_dc(idx_fix(idx_temp)),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
    hold on
    plot(time(idx_fix(idx_temp))-t1,Y_clean_dc(idx_fix(idx_temp)),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
end

% find all saccades
%inds = find(edf.events.sac_dc.trial == trial_property);
inds = find(edf.events.sac_dc.trial == trial_property & (edf.events.sac_dc.msg_srt == 11 | edf.events.sac_dc.msg_srt == 12));
sac_on = edf.events.sac_dc.ind_srt(inds);
sac_off = edf.events.sac_dc.ind_end(inds);
if strcmp(sac_property,'On')
    for jj = 1:length(inds)
        ind_sac = sac_on(jj):sac_off(jj);
        plot(time(ind_sac)-t1,X_clean_dc(ind_sac),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
        plot(time(ind_sac)-t1,Y_clean_dc(ind_sac),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
    end
end

% plot fixations
% find all fixations
inds_fix = find(edf.events.fix.trial == trial_property);
fix_on = edf.events.fix.ind_srt(inds_fix);
fix_off = edf.events.fix.ind_end(inds_fix);
if strcmp(fix_property,'On')
    for jj = 1:length(inds_fix)
        ind_fix = fix_on(jj):fix_off(jj);
        plot(time(ind_fix)-t1,X_clean_dc(ind_fix),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
        plot(time(ind_fix)-t1,Y_clean_dc(ind_fix),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
    end
end
box off
%xlim([1 17.5])
%ylim([-7 5])
xlim([1 16])
ylim([-3 5])
legend([p1,p2,p3,p4,p5,p6,p7,p8],{'Target X','Target Y', ...
    'Gaze X (raw)','Gaze Y (raw)', ...
    'Gaze X (cleaned)','Gaze Y (cleaned)' ...
    'Fixation X','Fixation Y'},'location','northeastoutside','box','off')
ax = gca;
ax.FontSize = 18;
print(f,'Fig1_eye_path_updated.png','-dpng','-r300')

%% Fig. 1E All trials

clear
close all
clc

% pick a random participant/run
subj = 'p16_1'; % p18_1

% plot
f = figure('Renderer', 'painters', 'Position', [10 10 1500 400])

for rr = 1:8
data_dir = fullfile('./v3_short',['sMGS' num2str(rr)],subj,'result');
file = dir(fullfile(data_dir,'*step6.mat'));
load(fullfile(file.folder,file.name),'edf','set')
edf_all{rr} = edf;

time = edf.samples.time;%(edf.samples.time - edf.samples.time(1))/1000; % time
X_clean_dc = edf.samples.x_deg_clean_dc(:,set.eye); % Gaze X
Y_clean_dc = edf.samples.y_deg_clean_dc(:,set.eye); % Gaze Y

for trial_property = 1:16
% index of a trial
ind_trial = find(edf.samples.trial == trial_property);
if edf.cal.resp_type(trial_property) < 6 | edf.cal.resp_type(trial_property) > 8
    continue;
end

% % start time of a trial
% srt_trial = time(ind_trial(1));
% srt_property = 0; % start time of a trial
% end_property = time(ind_trial(end))-time(ind_trial(1)); % end time of a trial
% 
% % plot specific time window within a trial
% t1 = srt_property + srt_trial; % srt_property +
% t2 = end_property + srt_trial; % end_property +

t1 = edf.events.msg.time(trial_property,2);
t2 = edf.events.msg.time(trial_property,end);

% find index of those times
ind = find(time >= t1 & time <= t2);

% find index of all the messages
n_msg = max(edf.samples.msg);
msgtime = edf.events.msg.time;
msgind = edf.events.msg.ind_srt;
msgtime_relative = msgtime - msgtime(:,2);

% delay offset
delay_on = edf.events.msg.time(trial_property,10);
delay_off = edf.events.msg.time(trial_property,11);
delay_dur = delay_off - delay_on;

% gaze x & y
% plot
x = X_clean_dc(ind);
y = Y_clean_dc(ind);

% if delay is short, let's add some nans to the end
delay_off_ind = msgind(trial_property,11) - msgind(trial_property,2);

% add nan
if round(delay_dur/1000) == 8
x = [x(1:delay_off_ind);nan(2000/(1000/edf.record.sample_rate),1);y(delay_off_ind:end)];
y = [y(1:delay_off_ind);nan(2000/(1000/edf.record.sample_rate),1);y(delay_off_ind:end)];
end

% % flip x if x is negative
% % flip y if y is positive
% time_resp = msgtime(trial_property,11):msgtime(trial_property,13);
% %time_resp = time_resp - msgtime(trial_property,1);
% ind_resp = find(time>time_resp(1) & time< time_resp(end)) - ind(1);
% x_mean_resp = mean(x(ind_resp),'omitnan');
% y_mean_resp = mean(y(ind_resp),'omitnan');
% if x_mean_resp < 0 % x is all positive
%     x = -x;
% end
% if y_mean_resp > 0 % y is all negative
%     y = -y;
% end
p1 = plot(1:length(x),x,'LineStyle',"-",'LineWidth',1.5,'Color',[0, 0.4470, 0.7410,1]);
hold on
p2 = plot(1:length(y),y,'LineStyle',"-",'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980,1]);
hold on

end
end

ax = gca;
ax.XTick = [0:2000:16000]/(1000/edf.record.sample_rate);
ax.XTickLabel = 0:2:16;

% plot messages
minY = -10;
maxY = 10;

line([1500, 1500]/(1000/edf.record.sample_rate), [minY,maxY], ...
    'LineStyle','--','Color','k','LineWidth',2);
hold on
line([5000, 5000]/(1000/edf.record.sample_rate), [minY,maxY], ...
    'LineStyle','--','Color','k','LineWidth',2);
hold on
line([15000, 15000]/(1000/edf.record.sample_rate), [minY,maxY], ...
    'LineStyle','--','Color','k','LineWidth',2);
hold on

box off
%xlim([1 17.5])
%ylim([-7 5])
xlim([0 16500]/(1000/edf.record.sample_rate))
ylim([minY maxY])
xlabel('Time (s)')
ylabel('Eye Position (dva)')
legend('off');
ax.FontSize = 18;

legend([p1,p2],{'X','Y'}, ...
    'location','northeastoutside','box','off')

print(f,'Fig1_eye_path_all_p16.png','-dpng','-r300')

