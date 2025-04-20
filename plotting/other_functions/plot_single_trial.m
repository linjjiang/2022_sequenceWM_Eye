% Plot preprocessing results (individual trials)
% By: Linjing Jiang
% Date: 04/08/2024

% Plot raw, cleaned, and drift corrected gaze locations
% Plot raw and cleaned pupil size time course


%%
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
% example, p18, sMGS1, trial 1, has good blink removal performance
%
output = '/mnt/Data1/eye_data/vswmda_eye/analysis/exp1';

%%

% check velocity distribution during fixation and delay period
% vel = edf.samples.vel(edf.samples.msg >= 2 & edf.samples.msg <= 3);
% figure;
% histogram(vel(vel<=500))

idx = find(edf.samples.trial == trial_property & (edf.samples.msg == 11 | edf.samples.msg == 12));
vel = edf.samples.vel_deg_clean_dc(idx,set.eye);
acc = edf.samples.acc_deg_clean_dc(idx,set.eye);
max(vel)
max(acc)
%%
close all
trial_property = 13; % trial
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

srt_property = 0;
end_property = time(ind_trial(end))-time(ind_trial(1));

% plot specific time window within a trial
t1 = srt_property + srt_trial; % srt_property +
t2 = end_property + srt_trial; % end_property +

% find index of those times
ind = find(time >= t1 & time <= t2);

% find index of all the messages
n_msg = max(edf.samples.msg);
%             for ii = 1:n_msg
%                 indm(ii,:) = find(edf.samples.msg == ii)';
%             end
indm = edf.events.msg.ind_srt;


% plot
f = figure('Renderer', 'painters', 'Position', [10 10 1500 800])
subplot(2,1,1)
plot(time(ind)-t1,p(ind),'LineStyle',"-","LineWidth",2,'color',[0,0.447,0.741,0.2]);
hold on
plot(time(ind)-t1,p_clean(ind),'LineStyle',"-","LineWidth",2,'color',[0,0.447,0.741,1]);
xlim([0 t2-t1])
ylim([min(p(ind)),max(p(ind))])
xlabel('Time (s)')
ylabel('Pupil Size')
title(['Pupil time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
legend('off');
%set(ZoomHandle,'Motion','horizontal')
% plot messages
hold on
if strcmp(msg_property,'On')
    for ii = trial_property %1:size(indm,1)
        for jj = 1:size(indm,2)
            h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[min(p(ind)),max(p(ind))],'LineStyle','--','Color','k');%cmap(ii,:)
            text(time(indm(ii,jj))-t1,max(p(ind))*0.9,num2str(edf.events.msg.txt(ii,jj))); % set.msg{jj}
        end
    end
end
box off
% gaze x & y
subplot(2,1,2)

% plot target locations

line([0 t2-t1],[stimX(taridx) stimX(taridx)],'LineStyle','--','Color',[0, 0.4470, 0.7410,1])
hold on
line([0 t2-t1],[stimY(taridx) stimY(taridx)],'LineStyle','--','Color',[0.8500, 0.3250, 0.0980,1])
hold on
line([0 t2-t1],[tarX tarX],'LineStyle','--','Color',[0, 0.4470, 0.7410,1])
hold on
line([0 t2-t1],[tarY tarY],'LineStyle','--','Color',[0.8500, 0.3250, 0.0980,1])
hold on


% minY and maxY
minY = min([X(ind);Y(ind);tarX;tarY;X_clean_dc(ind);Y_clean_dc(ind)])-0.5;
maxY = max([X(ind);Y(ind);tarX;tarY;X_clean_dc(ind);Y_clean_dc(ind)])+0.5;
% gaze locations
plot(time(ind)-t1,X(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,0.2]);
hold on
plot(time(ind)-t1,Y(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,0.2]);
% plot(time(ind)-t1,X_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
% plot(time(ind)-t1,Y_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
plot(time(ind)-t1,X_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,1]);
plot(time(ind)-t1,Y_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,1]);
xlim([0 t2-t1])
ylim([minY maxY])
xlabel('Time (s)')
ylabel('Eye Position (deg)')
legend('off');
title(['Gaze time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
hold on

% plot messages
if strcmp(msg_property,'On')
    for ii = trial_property%1:size(indm,1)
        for jj = [3,10,11]%1:size(indm,2)
            h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[minY,maxY],'LineStyle','--','Color','k');% cmap(ii,:)
           % text(time(indm(ii,jj))-t1,maxY*0.9,0,num2str(edf.events.msg.txt(ii,jj)));
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
    plot(time(pfix_on:pfix_off)-t1,X_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
    plot(time(pfix_on:pfix_off)-t1,Y_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
end




% plot fixations used to perform drift correction
idx_fix = find(edf.samples.baseline_dc == 1 & ... % is a baseline
    edf.samples.trial == trial_property); 

% idx_fix1 = find(edf.samples.trial == trial_property & ... % trial number equals ii
%                edf.samples.msg == 2 & ... % during fixation period (message 2-3)
%               abs(edf.samples.x_deg_clean(:,set.eye)) <= set.noise.dc_threshold & ... % x position less than 3 dva
%               abs(edf.samples.y_deg_clean(:,set.eye)) <= set.noise.dc_threshold); % y position less than 3 dva

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
xlim([1 17.5])
ylim([-8 4])

%saveas(f,fullfile(output,[filename,'_',run,'_trial_',num2str(trial_property),'_primary_any.jpg']))
% 
%% Plot most accurate saccade
% plot
f = figure('Renderer', 'painters', 'Position', [10 10 1500 800])
subplot(2,1,1)
plot(time(ind)-t1,p(ind),'LineStyle',"-","LineWidth",2,'color',[0,0.447,0.741,0.2]);
hold on
plot(time(ind)-t1,p_clean(ind),'LineStyle',"-","LineWidth",2,'color',[0,0.447,0.741,1]);
xlim([0 t2-t1])
ylim([min(p(ind)),max(p(ind))])
xlabel('Time (s)')
ylabel('Pupil Size')
title(['Pupil time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
legend('off');
%set(ZoomHandle,'Motion','horizontal')
% plot messages
hold on
if strcmp(msg_property,'On')
    for ii = trial_property %1:size(indm,1)
        for jj = 1:size(indm,2)
            h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[min(p(ind)),max(p(ind))],'LineStyle','--','Color','k');%cmap(ii,:)
            text(time(indm(ii,jj))-t1,max(p(ind))*0.9,num2str(edf.events.msg.txt(ii,jj))); % set.msg{jj}
        end
    end
end
box off
% gaze x & y
subplot(2,1,2)

% plot target locations

line([0 t2-t1],[tarX tarX],'LineStyle','--','Color',[0, 0.4470, 0.7410,1])
hold on
line([0 t2-t1],[tarY tarY],'LineStyle','--','Color',[0.8500, 0.3250, 0.0980,1])
hold on

% minY and maxY
minY = min([X(ind);Y(ind);tarX;tarY;X_clean_dc(ind);Y_clean_dc(ind)])-0.5;
maxY = max([X(ind);Y(ind);tarX;tarY;X_clean_dc(ind);Y_clean_dc(ind)])+0.5;
% gaze locations
plot(time(ind)-t1,X(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,0.2]);
hold on
plot(time(ind)-t1,Y(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,0.2]);
% plot(time(ind)-t1,X_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
% plot(time(ind)-t1,Y_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
plot(time(ind)-t1,X_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,1]);
plot(time(ind)-t1,Y_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,1]);
xlim([0 t2-t1])
ylim([minY maxY])
xlabel('Time (s)')
ylabel('Eye Position (deg)')
legend('off');
title(['Gaze time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
hold on

% plot messages
if strcmp(msg_property,'On')
    for ii = trial_property%1:size(indm,1)
        for jj = 1:size(indm,2)
            h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[minY,maxY],'LineStyle','--','Color','k');% cmap(ii,:)
            text(time(indm(ii,jj))-t1,maxY*0.9,0,num2str(edf.events.msg.txt(ii,jj)));
        end
    end
end
hold on


% plot saccades

% plot most accurate saccade

ind_psac = edf.cal.acc_sac_tar_ind(trial_property);
if ~isnan(ind_psac)
    psac_on = edf.events.sac_dc.ind_srt(ind_psac);
    psac_off = edf.events.sac_dc.ind_end(ind_psac);
    plot(time(psac_on:psac_off)-t1,X_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
    plot(time(psac_on:psac_off)-t1,Y_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);

% endpoint fixations
    pfix_on = edf.events.fix_dc.ind_srt(edf.events.sac_dc.endfix(ind_psac));
    pfix_off = edf.events.fix_dc.ind_end(edf.events.sac_dc.endfix(ind_psac));
    plot(time(pfix_on:pfix_off)-t1,X_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
    plot(time(pfix_on:pfix_off)-t1,Y_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);

end

% plot fixations used to perform drift correction
idx_fix = find(edf.samples.baseline_dc == 1 & ... % is a baseline
    edf.samples.trial == trial_property); 

% idx_fix1 = find(edf.samples.trial == trial_property & ... % trial number equals ii
%                edf.samples.msg == 2 & ... % during fixation period (message 2-3)
%               abs(edf.samples.x_deg_clean(:,set.eye)) <= set.noise.dc_threshold & ... % x position less than 3 dva
%               abs(edf.samples.y_deg_clean(:,set.eye)) <= set.noise.dc_threshold); % y position less than 3 dva

dif_fix = diff(idx_fix);
idx_fix_sep = find(dif_fix > 1); idx_fix_sep = [0;idx_fix_sep;length(idx_fix)];
for xx = 1:(length(idx_fix_sep)-1)
    idx_temp = (idx_fix_sep(xx)+1):(idx_fix_sep(xx+1));
    plot(time(idx_fix(idx_temp))-t1,X_clean_dc(idx_fix(idx_temp)),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
    hold on
    plot(time(idx_fix(idx_temp))-t1,Y_clean_dc(idx_fix(idx_temp)),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
end

%     saveas(f,fullfile(output,[filename,'_',run,'_trial_',num2str(trial_property),'_acc_tar.jpg']))
% 
% %% plot the primary saccade that enters the target AOI
% % plot
% f = figure('Renderer', 'painters', 'Position', [10 10 1500 800])
% subplot(2,1,1)
% plot(time(ind)-t1,p(ind),'LineStyle',"-","LineWidth",2,'color',[0,0.447,0.741,0.2]);
% hold on
% plot(time(ind)-t1,p_clean(ind),'LineStyle',"-","LineWidth",2,'color',[0,0.447,0.741,1]);
% xlim([0 t2-t1])
% ylim([min(p(ind)),max(p(ind))])
% xlabel('Time (s)')
% ylabel('Pupil Size')
% title(['Pupil time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
% legend('off');
% %set(ZoomHandle,'Motion','horizontal')
% % plot messages
% hold on
% if strcmp(msg_property,'On')
%     for ii = trial_property %1:size(indm,1)
%         for jj = 1:size(indm,2)
%             h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[min(p(ind)),max(p(ind))],'LineStyle','--','Color','k');%cmap(ii,:)
%             text(time(indm(ii,jj))-t1,max(p(ind))*0.9,num2str(edf.events.msg.txt(ii,jj))); % set.msg{jj}
%         end
%     end
% end
% box off
% % gaze x & y
% subplot(2,1,2)
% 
% % plot target locations
% 
% line([0 t2-t1],[tarX tarX],'LineStyle','--','Color',[0, 0.4470, 0.7410,1])
% hold on
% line([0 t2-t1],[tarY tarY],'LineStyle','--','Color',[0.8500, 0.3250, 0.0980,1])
% hold on
% 
% % minY and maxY
% minY = min([X(ind);Y(ind);tarX;tarY])-0.5;
% maxY = max([X(ind);Y(ind);tarX;tarY])+0.5;
% % gaze locations
% plot(time(ind)-t1,X(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,0.2]);
% hold on
% plot(time(ind)-t1,Y(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,0.2]);
% % plot(time(ind)-t1,X_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
% % plot(time(ind)-t1,Y_clean(ind),'LineStyle',"--",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
% plot(time(ind)-t1,X_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0, 0.4470, 0.7410,1]);
% plot(time(ind)-t1,Y_clean_dc(ind),'LineStyle',"-",'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980,1]);
% xlim([0 t2-t1])
% ylim([minY maxY])
% xlabel('Time (s)')
% ylabel('Eye Position (deg)')
% legend('off');
% title(['Gaze time course: trial ',num2str(trial_property),' ',erase(filename,'_1'),' ',run])
% hold on
% 
% % plot messages
% if strcmp(msg_property,'On')
%     for ii = trial_property%1:size(indm,1)
%         for jj = 1:size(indm,2)
%             h = line([time(indm(ii,jj))-t1,time(indm(ii,jj))-t1],[minY,maxY],'LineStyle','--','Color','k');% cmap(ii,:)
%             text(time(indm(ii,jj))-t1,maxY*0.9,0,num2str(edf.events.msg.txt(ii,jj)));
%         end
%     end
% end
% hold on
% 
% 
% % plot saccades
% 
% % plot the primary saccade that enters the target AOI
% ind_psac = edf.cal.primary_sac_tar_ind(trial_property);
% if ~isnan(ind_psac)
%     psac_on = edf.events.sac_dc.ind_srt(ind_psac);
%     psac_off = edf.events.sac_dc.ind_end(ind_psac);
%     plot(time(psac_on:psac_off)-t1,X_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
%     plot(time(psac_on:psac_off)-t1,Y_clean_dc(psac_on:psac_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
% % endpoint fixations
%     pfix_on = edf.events.fix_dc.ind_srt(edf.events.sac_dc.endfix(ind_psac));
%     pfix_off = edf.events.fix_dc.ind_end(edf.events.sac_dc.endfix(ind_psac));
%     plot(time(pfix_on:pfix_off)-t1,X_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
%     plot(time(pfix_on:pfix_off)-t1,Y_clean_dc(pfix_on:pfix_off),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
% 
% end
% 
% 
% % plot fixations used to perform drift correction
% idx_fix = find(edf.samples.baseline_dc == 1 & ... % is a baseline
%     edf.samples.trial == trial_property); 
% 
% % idx_fix1 = find(edf.samples.trial == trial_property & ... % trial number equals ii
% %                edf.samples.msg == 2 & ... % during fixation period (message 2-3)
% %               abs(edf.samples.x_deg_clean(:,set.eye)) <= set.noise.dc_threshold & ... % x position less than 3 dva
% %               abs(edf.samples.y_deg_clean(:,set.eye)) <= set.noise.dc_threshold); % y position less than 3 dva
% 
% dif_fix = diff(idx_fix);
% idx_fix_sep = find(dif_fix > 1); idx_fix_sep = [0;idx_fix_sep;length(idx_fix)];
% for xx = 1:(length(idx_fix_sep)-1)
%     idx_temp = (idx_fix_sep(xx)+1):(idx_fix_sep(xx+1));
%     plot(time(idx_fix(idx_temp))-t1,X_clean_dc(idx_fix(idx_temp)),'LineWidth',6,'Color',[0, 0.4470, 0.7410]);
%     hold on
%     plot(time(idx_fix(idx_temp))-t1,Y_clean_dc(idx_fix(idx_temp)),'LineWidth',6,'Color',[0.8500, 0.3250, 0.0980]);
% end%    saveas(f,fullfile(output,[filename,'_',run,'_trial_',num2str(trial_property),'_primary_tar.jpg']))
