%%Alternative Fig 3 Exp 1 
% Plot modeling results

% By: Linjing Jiang
% Date: 04/22/2024

%% Load the data
clear
close all
clc

% output directory
output = './analysis/exp2';


cd(output)

% all the csv files
% p_err_wide.csv: primary saccade error (for all stimuli)
% p_rt_wide.csv: primary saccade latency (for all stimuli)
% pt_err_wide.csv: primary saccade error (for target)
% pt_rt_wide.csv: primary saccade latency (for target)
% a_err_wide.csv: most accurate saccade error (for target)

f = figure('Renderer', 'painters', 'Position', [10 10 800 1500]) %1350
fz = 14; 
markersz = 8;
linewz = 1.5;
tiledlayout(3,2)


%% Plot sd
nexttile

tbl = readtable(fullfile(output,'p_sd_out.csv'));
cue = tbl.cue;
data = [tbl{cue==0,1:4} [tbl{cue==1,1:4};nan(5,4)]];

m = mean(data,'omitnan');
sd = [std(data(:,1:4),0,1,'omitnan')/sqrt(sum(cue==0)) ...
    std(data(:,5:8),0,1,'omitnan')/sqrt(sum(cue==1))];

% effects of serial position and cue
% the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% respectively
% scatterplot
xCenter = ones(size(data)).*[1:4,6:9];
xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
xCenter(isnan(data)) = NaN;
x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
%colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
colors = repmat([0 0 0],size(data,2)*size(data,1),1);
%scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% for ii = 1:length(x1)
% scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% hold on
% end
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',linewz,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',markersz) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',linewz,'Color',[0 0 0],'Marker','o','MarkerSize',markersz) % [217,95,2]./255
xlim([0 5])
ylim([1.2 2])
% legend('Quadrant cue','Order cue','box','off')
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''}, ...
    'YTick',[1.2:0.2:2],'YTickLabel',{'1.2','1.4','1.6','1.8','2.0'})
ylabel('Memory imprecision (dva)')
%xlabel('Serial position')
%title('First saccade')
box off


%% Swap rate
nexttile
tbl = readtable(fullfile(output,'p_beta_out.csv'));
cue = tbl.cue;
data = [tbl{cue==0,1:4} [tbl{cue==1,1:4};nan(5,4)]];

m = mean(data,'omitnan');
sd = [std(data(:,1:4),0,1,'omitnan')/sqrt(sum(cue==0)) ...
    std(data(:,5:8),0,1,'omitnan')/sqrt(sum(cue==1))];

% effects of serial position and cue
% the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% respectively
% scatterplot
xCenter = ones(size(data)).*[1:4,6:9];
xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
xCenter(isnan(data)) = NaN;
x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
%colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
colors = repmat([0 0 0],size(data,2)*size(data,1),1);
%scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% for ii = 1:length(x1)
% scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% hold on
% end
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',linewz,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',markersz) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',linewz,'Color',[0 0 0],'Marker','o','MarkerSize',markersz) % [217,95,2]./255
xlim([0 5])
%ylim([1.3 3.2])
legend('Quadrant cue','Order cue','box','off','FontSize',fz) %(N = 65)
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''}, ...
    'YTick',[0:0.1:0.2],'YTickLabel',{'0','0.1','0.2'})
ylabel('Swap rate')
%xlabel('Serial position')
%title('Swap rate')
box off


%% Target rate
nexttile

tbl = readtable(fullfile(output,'p_alpha_out.csv'));
cue = tbl.cue;
data = [tbl{cue==0,1:4} [tbl{cue==1,1:4};nan(5,4)]];

m = mean(data,'omitnan');
sd = [std(data(:,1:4),0,1,'omitnan')/sqrt(sum(cue==0)) ...
    std(data(:,5:8),0,1,'omitnan')/sqrt(sum(cue==1))];

% effects of serial position and cue
% the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% respectively
% scatterplot
xCenter = ones(size(data)).*[1:4,6:9];
xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
xCenter(isnan(data)) = NaN;
x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
%colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
colors = repmat([0 0 0],size(data,2)*size(data,1),1);
%scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% for ii = 1:length(x1)
% scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% hold on
% end
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',linewz,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',markersz) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',linewz,'Color',[0 0 0],'Marker','o','MarkerSize',markersz) % [217,95,2]./255
xlim([0 5])
%ylim([550 860])
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''}, ...
    'YTick',[0.8:0.1:1],'YTickLabel',{'0.8','0.9','1.0'})
ylabel('Target rate')
%xlabel('Serial position')
box off

%% Guessing rate
nexttile
tbl = readtable(fullfile(output,'p_gamma_out.csv'));
cue = tbl.cue;
data = [tbl{cue==0,1:4} [tbl{cue==1,1:4};nan(5,4)]];

m = mean(data,'omitnan');
sd = [std(data(:,1:4),0,1,'omitnan')/sqrt(sum(cue==0)) ...
    std(data(:,5:8),0,1,'omitnan')/sqrt(sum(cue==1))];

% effects of serial position and cue
% the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% respectively
% scatterplot
%xCenter = ones(size(data)).*[1:4,6:9];
%xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
%xCenter(isnan(data)) = NaN;
%x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
%%colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
%colors = repmat([0 0 0],size(data,2)*size(data,1),1);
%scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% for ii = 1:length(x1)
% scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% hold on
% end
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',linewz,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',markersz) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',linewz,'Color',[0 0 0],'Marker','o','MarkerSize',markersz) % [217,95,2]./255
xlim([0 5])
ylim([0 0.06])
%legend('3-dva','5-dva','box','off')
% set(gca,'FontSize',20,'XTick',1:9,'XTickLabel',{'1st','2nd','3rd','4th','','1st','2nd','3rd','4th'})
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''},...
    'YTick',0:0.02:0.06)
ylabel('Guessing rate')
%xlabel('Serial position')
%title('First saccade to target')
box off

%% Responded serial position
cd(output)
load('stats_exp2_quad.mat') % stats1
load('stats_exp2_ord.mat') % stats2

nexttile
final_sind1 = stats1.final_sind;
final_sind2 = stats2.final_sind;

% transposition error % for all the not nan trials
p_ord_err = [stats1.p_ord_err_ord(final_sind1,:,:) [stats2.p_ord_err_ord(final_sind2,:,:);nan(5,64,8)]];

% the responded serial position
for pp = 1:4
p_sp(:,:,[pp,pp+4]) = p_ord_err(:,:,[pp,pp+4])+pp; 
end

% plot
% first participant, first serial position
p_sp_perc = nan(size(p_sp,1),4,8);
for ss = 1:length(final_sind1) % for each participant
    for pp = 1:4 % for each serial position
temp = squeeze(p_sp(ss,:,pp)); % get 16 trials
tr_perc = [sum(temp==1) sum(temp==2) ...
    sum(temp==3) sum(temp==4)]; %/16*100 %sum(~isnan(temp))
p_sp_perc(ss,:,pp) = tr_perc;
    end
end

for ss = 1:length(final_sind2) % for each participant
    for pp = 5:8 % for each serial position
temp = squeeze(p_sp(ss,:,pp)); % get 16 trials
tr_perc = [sum(temp==1) sum(temp==2) ...
    sum(temp==3) sum(temp==4)]; %/16*100 %sum(~isnan(temp))
p_sp_perc(ss,:,pp) = tr_perc;
    end
end

% plot the percentage across participants
%markers = {'o','x','+','*'};
%lstyle = {':','-.','--','-'};
lstyle = {'-','-','-','-'};
for pp = 1:4
% for ss = 1:size(p_ord_err,1)
% plot([1:4],p_sp_perc(ss,:,pp),'color',[0.5 0.5 0.5],'LineWidth',1) %[27,158,119,50]./255
% hold on
% plot([1:4],p_sp_perc(ss,:,pp+4),'color',[0 0 0],'LineWidth',1) %[217,95,2,50]./255
% hold on
% end
errorbar([1:4],mean(p_sp_perc(:,:,pp),1,'omitnan'),std(p_sp_perc(:,:,pp),0,1,'omitnan')/sqrt(length(final_sind1)),'color',[0.5 0.5 0.5],'LineWidth',linewz,'LineStyle','--','Marker','s','MarkerSize',markersz)
hold on
errorbar([1:4],mean(p_sp_perc(:,:,pp+4),1,'omitnan'),std(p_sp_perc(:,:,pp+4),0,1,'omitnan')/sqrt(length(final_sind2)),'color',[0 0 0],'LineWidth',linewz,'LineStyle','-','Marker','o','MarkerSize',markersz)
hold on
end
box off
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
xlim([0 5])
xlabel('Serial position')
ylabel('Number of trials')
box off

%% Transposition error - cue by serial position
nexttile
tbl = readtable(fullfile(output,'p_orderr_out.csv'));
cue = tbl.cue;
data = [tbl{cue==0,1:4} [tbl{cue==1,1:4};nan(5,4)]];

m = mean(data,'omitnan');
sd = [std(data(:,1:4),0,1,'omitnan')/sqrt(sum(cue==0)) ...
    std(data(:,5:8),0,1,'omitnan')/sqrt(sum(cue==1))];

% effects of serial position and cue
% the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% respectively
% scatterplot
xCenter = ones(size(data)).*[1:4,6:9];
xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
xCenter(isnan(data)) = NaN;
x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
%colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
colors = repmat([0 0 0],size(data,2)*size(data,1),1);
%scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% for ii = 1:length(x1)
% scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% hold on
% end
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',linewz,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',markersz) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',linewz,'Color',[0 0 0],'Marker','o','MarkerSize',markersz) % [217,95,2]./255
xlim([0 5])
%legend('Quadrant cue','Order cue','box','off') %'location','northeastoutside'
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
ylabel('Transposition error (%)')
xlabel('Serial position')
box off

print(f,'Fig5_exp2_updated.png','-dpng','-r300')

%saveas(f,fullfile(output,['Fig5_exp2.jpg']))
%print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['Fig5_exp2.eps']))

% %% Plot mu
% 
% tbl = readtable(fullfile(output,'p_mu_wide.csv'));
% data = tbl{:,1:8};
% 
% m = mean(data,'omitnan');
% sd = std(data,0,1,'omitnan')/sqrt(size(data,1));
% 
% % effects of serial position and cue
% % the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% % respectively
% f = figure('Renderer', 'painters', 'Position', [10 10 800 600])
% % scatterplot
% xCenter = ones(size(data)).*[1:4,6:9];
% xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
% xCenter(isnan(data)) = NaN;
% x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
% %colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
% colors = repmat([0 0 0],size(data,2)*size(data,1),1);
% %scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% % for ii = 1:length(x1)
% % scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% % hold on
% % end
% errorbar(1:4,m(1:4),sd(1:4),'LineWidth',4,'Color',[27,158,119]./255,'Marker','^','MarkerSize',10) %
% hold on
% errorbar(1:4,m(5:8),sd(5:8),'LineWidth',4,'Color',[217,95,2]./255,'Marker','o','MarkerSize',10) %
% xlim([0 10])
% %legend('3-dva','5-dva','box','off')
% set(gca,'FontSize',20,'XTick',1:9,'XTickLabel',{'1st','2nd','3rd','4th','','1st','2nd','3rd','4th'})
% ylabel('Modeled systematic bias (pixel)')
% xlabel('Serial position')
% box off
% saveas(f,fullfile(output,['mu.jpg']))
% print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['mu.eps']))
% 
% %% Plot sd
% 
% tbl = readtable(fullfile(output,'p_sd_wide.csv'));
% data = tbl{:,1:8};
% 
% m = mean(data,'omitnan');
% sd = std(data,0,1,'omitnan')/sqrt(size(data,1));
% 
% % effects of serial position and cue
% % the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% % respectively
% f = figure('Renderer', 'painters', 'Position', [10 10 800 600])
% % scatterplot
% xCenter = ones(size(data)).*[1:4,6:9];
% xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
% xCenter(isnan(data)) = NaN;
% x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
% %colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
% colors = repmat([0 0 0],size(data,2)*size(data,1),1);
% %scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% % for ii = 1:length(x1)
% % scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% % hold on
% % end
% errorbar(1:4,m(1:4),sd(1:4),'LineWidth',4,'Color',[27,158,119]./255,'Marker','^','MarkerSize',10) %
% hold on
% errorbar(1:4,m(5:8),sd(5:8),'LineWidth',4,'Color',[217,95,2]./255,'Marker','o','MarkerSize',10) %
% xlim([0 10])
% %legend('3-dva','5-dva','box','off')
% set(gca,'FontSize',20,'XTick',1:9,'XTickLabel',{'1st','2nd','3rd','4th','','1st','2nd','3rd','4th'})
% ylabel('Modeled memory imprecision (pixel)')
% xlabel('Serial position')
% box off
% saveas(f,fullfile(output,['sd.jpg']))
% print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['sd.eps']))
% 
% %% Plot beta
% 
% tbl = readtable(fullfile(output,'p_beta_wide.csv'));
% data = tbl{:,1:8};
% 
% m = mean(data,'omitnan');
% sd = std(data,0,1,'omitnan')/sqrt(size(data,1));
% 
% % effects of serial position and cue
% % the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% % respectively
% f = figure('Renderer', 'painters', 'Position', [10 10 800 600])
% % scatterplot
% xCenter = ones(size(data)).*[1:4,6:9];
% xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
% xCenter(isnan(data)) = NaN;
% x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
% %colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
% colors = repmat([0 0 0],size(data,2)*size(data,1),1);
% %scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% % for ii = 1:length(x1)
% % scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% % hold on
% % end
% errorbar(1:4,m(1:4),sd(1:4),'LineWidth',4,'Color',[27,158,119]./255,'Marker','^','MarkerSize',10) %
% hold on
% errorbar(1:4,m(5:8),sd(5:8),'LineWidth',4,'Color',[217,95,2]./255,'Marker','o','MarkerSize',10) %
% xlim([0 10])
% %legend('3-dva','5-dva','box','off')
% set(gca,'FontSize',20,'XTick',1:9,'XTickLabel',{'1st','2nd','3rd','4th','','1st','2nd','3rd','4th'})
% ylabel('Modeled swap rate')
% xlabel('Serial position')
% box off
% saveas(f,fullfile(output,['beta.jpg']))
% print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['beta.eps']))
% 
% %% Plot gamma
% tbl = readtable(fullfile(output,'p_gamma_wide.csv'));
% data = tbl{:,1:8};
% 
% m = mean(data,'omitnan');
% sd = std(data,0,1,'omitnan')/sqrt(size(data,1));
% 
% % effects of serial position and cue
% % the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% % respectively
% f = figure('Renderer', 'painters', 'Position', [10 10 800 600])
% % scatterplot
% xCenter = ones(size(data)).*[1:4,6:9];
% xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
% xCenter(isnan(data)) = NaN;
% x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
% %colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
% colors = repmat([0 0 0],size(data,2)*size(data,1),1);
% %scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% % for ii = 1:length(x1)
% % scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% % hold on
% % end
% errorbar(1:4,m(1:4),sd(1:4),'LineWidth',4,'Color',[27,158,119]./255,'Marker','^','MarkerSize',10) %
% hold on
% errorbar(1:4,m(5:8),sd(5:8),'LineWidth',4,'Color',[217,95,2]./255,'Marker','o','MarkerSize',10) %
% xlim([0 10])
% %legend('3-dva','5-dva','box','off')
% set(gca,'FontSize',20,'XTick',1:9,'XTickLabel',{'1st','2nd','3rd','4th','','1st','2nd','3rd','4th'})
% ylabel('Modeled guessing rate')
% xlabel('Serial position')
% box off
% saveas(f,fullfile(output,['gamma.jpg']))
% print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['gamma.eps']))

% %% Transposition error
% cd(output)
% load('stats_exp1_longdelay.mat') % stats1
% load('stats_exp1_shortdelay.mat') % stats2
% final_sind1 = stats1.final_sind;
% final_sind2 = stats2.final_sind;
% 
% % transposition error % for all the not nan trials
% p_ord_err = [stats1.p_ord_err_ord(final_sind1,:,:);stats2.p_ord_err_ord(final_sind2,:,:)];
% 
% % the responded serial position
% for pp = 1:4
% p_sp(:,:,[pp,pp+4]) = p_ord_err(:,:,[pp,pp+4])+pp; 
% end
% 
% % plot
% % first participant, first serial position
% for ss = 1:size(p_ord_err,1)
%     for pp = 1:size(p_ord_err,3)
% temp = squeeze(p_sp(ss,:,pp));
% tr_perc = [sum(temp==1)/16*100 sum(temp==2)/16*100 ...
%     sum(temp==3)/16*100 sum(temp==4)/16*100]; %sum(~isnan(temp))
% p_sp_perc(ss,:,pp) = tr_perc;
%     end
% end
% 
% % plot the percentage across participants
% close all
% f = figure('Renderer', 'painters', 'Position', [10 10 800 600])
% %markers = {'o','x','+','*'};
% %lstyle = {':','-.','--','-'};
% lstyle = {'-','-','-','-'};
% for pp = 1:4
% for ss = 1:size(p_ord_err,1)
% plot([1:4],p_sp_perc(ss,:,pp),'color',[27,158,119,50]./255,'LineWidth',1)
% hold on
% plot([1:4],p_sp_perc(ss,:,pp+4),'color',[217,95,2,50]./255,'LineWidth',1)
% hold on
% end
% plot([1:4],mean(p_sp_perc(:,:,pp),1),'color',[27,158,119,255]./255,'LineWidth',4,'LineStyle',lstyle{pp},'Marker','^','MarkerSize',10)
% hold on
% plot([1:4],mean(p_sp_perc(:,:,pp+4),1),'color',[217,95,2,255]./255,'LineWidth',4,'LineStyle',lstyle{pp},'Marker','o','MarkerSize',10)
% hold on
% end
% box off
% set(gca,'FontSize',20,'XTick',1:4)
% xlabel('Responded serial position')
% ylabel('Percentage of trials (%)')
% box off
% saveas(f,fullfile(output,['perc_order_resp.jpg']))
% print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['perc_order_resp.eps']))
% 
% %% Transposition error - cue by serial position
% tbl = readtable(fullfile(output,'p_orderr_wide.csv'));
% data = tbl{:,1:8};
% 
% m = mean(data,'omitnan');
% sd = std(data,0,1,'omitnan')/sqrt(size(data,1));
% 
% % effects of serial position and cue
% % the figure has 8 columns, 1-4 serial positions for quadrant and order cue
% % respectively
% f = figure('Renderer', 'painters', 'Position', [10 10 500 600])
% % scatterplot
% xCenter = ones(size(data)).*[1:4,6:9];
% xCenter = xCenter + (-1 + 2* rand(size(data)))*0.1;
% xCenter(isnan(data)) = NaN;
% x1 = reshape(xCenter,1,[]); y1 = reshape(data,1,[]);
% %colors = repmat(distinguishable_colors(size(data,1),{'w','k'}),size(data,2),1);
% colors = repmat([0 0 0],size(data,2)*size(data,1),1);
% %scatter(x1,y1,'lineWidth',0.1,'markerFaceColor',[0.8 0.8 0.8],'markerEdgeColor',[0.2 0.2 0.2])
% % for ii = 1:length(x1)
% % scatter(x1(ii),y1(ii),'lineWidth',0.1,'MarkerFaceColor',colors(ii,:),'markerEdgeColor',colors(ii,:),'SizeData',5)
% % hold on
% % end
% errorbar(1:4,m(1:4),sd(1:4),'LineWidth',4,'Color',[27,158,119]./255,'Marker','^','MarkerSize',10) %
% hold on
% errorbar(1:4,m(5:8),sd(5:8),'LineWidth',4,'Color',[217,95,2]./255,'Marker','o','MarkerSize',10) %
% xlim([0 5])
% %legend('3-dva','5-dva','box','off')
% set(gca,'FontSize',20,'XTick',0:5,'XTickLabel',{'','1st','2nd','3rd','4th',''})
% ylabel('Transposition error (%)')
% xlabel('Serial position')
% box off
% saveas(f,fullfile(output,['trans_error_for_all.jpg']))
% print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['trans_error_for_all.eps']))


