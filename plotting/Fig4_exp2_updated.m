%% Plot the mean not box plot

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

f = figure('Renderer', 'painters', 'Position', [10 10 1100 1000])
fz = 20; 
tiledlayout(2,2)

%% Plot primary saccade errors (for all stimuli)
nexttile

tbl = readtable(fullfile(output,'p_err_out.csv'));
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
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',10) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',2,'Color',[0 0 0],'Marker','o','MarkerSize',10) % [217,95,2]./255
xlim([0 5])
ylim([1.5 3.5])
% legend('Quadrant cue','Order cue','box','off')
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''}, ...
    'YTick',[1.5:0.5:3.5],'YTickLabel',{'1.5','2.0','2.5','3.0','3.5'})
ylabel('Error (dva)')
%xlabel('Serial position')
%title('First saccade')
box off


%% Plot most accurate saccade errors (for target stimuli)
nexttile
tbl = readtable(fullfile(output,'a_err_out.csv'));
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
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',10) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',2,'Color',[0 0 0],'Marker','o','MarkerSize',10) % [217,95,2]./255
xlim([0 5])
ylim([1.2 2])
legend('Quadrant cue','Order cue','box','off','FontSize',fz) %(N = 65)
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''}, ...
    'YTick',[1.2:0.2:2],'YTickLabel',{'1.2','1.4','1.6','1.8','2.0'})
%ylabel('Error (dva)')
%xlabel('Serial position')
%title('Closest saccade to target')

box off

%% Plot primary saccade latency (for all stimuli)
nexttile

tbl = readtable(fullfile(output,'p_rt_out.csv'));
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
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',10) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',2,'Color',[0 0 0],'Marker','o','MarkerSize',10) % [217,95,2]./255
xlim([0 5])
ylim([450 950])
%legend('3-dva','5-dva','box','off')
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
ylabel('Latency (ms)')
xlabel('Serial position')
box off


%% Latency for most accurate saccade
nexttile
tbl = readtable(fullfile(output,'a_rt_out.csv'));
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
errorbar(1:4,m(1:4),sd(1:4),'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5],'Marker','s','MarkerSize',10) % [27,158,119]./255
hold on
errorbar(1:4,m(5:8),sd(5:8),'LineWidth',2,'Color',[0 0 0],'Marker','o','MarkerSize',10) % [217,95,2]./255
xlim([0 5])
ylim([450 950])
%legend('3-dva','5-dva','box','off')
set(gca,'FontSize',fz,'XTick',0:5,'XTickLabel',{'','1','2','3','4',''})
%ylabel('Latency (ms)')
xlabel('Serial position')
box off

print(f,'Fig4_exp2_updated.png','-dpng','-r300')
%print(f, '-depsc2', '-r300', '-tiff', '-loose', fullfile(output,['Fig4_exp2.eps']))

