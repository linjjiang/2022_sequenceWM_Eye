% Get demographics for Linjing's thesis - behavioral experiment
% VerDate: 07/25/2024
clear 
close all
clc

%% load subject id from experiment 1
cd('./analysis/exp1')
load('stats_exp1_longdelay.mat')
load('stats_exp1_shortdelay.mat')

sub_id_exp1 = [stats1.subs_id(stats1.final_sind) ...
    stats2.subs_id(stats2.final_sind)];
% In total 65 participants

%% load subject id from experiment 2
% cd('./analysis/exp2')
% load('stats_exp2_quad.mat')
% load('stats_exp2_ord.mat')
% 
% sub_id_exp2 = [stats1.subs_id(stats1.final_sind) ...
%     stats2.subs_id(stats2.final_sind)];

sub_id_exp2 = {'s108_1','s112_1','s116_1','s125_1','s127_1',...
's129_1','s135_1','s139_1','s145_1','s147_1','s149_1','s153_1',...
's155_1','s159_1','s161_1','s163_1','s107_1','s115_1','s122_1',...
's126_1','s128_1','s136_1','s142_1','s148_1','s150_1','s152_1',...
's164_1'};

%% Remove the _1 from the id
sub_id_exp1 = strrep(sub_id_exp1,'_1','');
sub_id_exp2 = strrep(sub_id_exp2,'_1','');


%% load demographics
cd('./analysis/demographics')

tbl1 = readtable(files(1).name,'VariableNamingRule','preserve');
tbl2 = readtable(files(2).name,'VariableNamingRule','preserve');

tbl2.Sex = cellfun(@(c) double(c~='F'),tbl2.Sex,'UniformOutput',false);
tbl2.Sex(cellfun(@isempty,tbl2.Sex)) = {NaN};

temp = cell2mat(tbl2.Sex);
tbl2.Sex = temp;

tbl = [tbl1;tbl2];

%% Calculate age and gender
age_exp1 = tbl.Age(ismember(tbl.ID,sub_id_exp1));
age_exp2 = tbl.Age(ismember(tbl.ID,sub_id_exp2));
gender_exp1 = tbl.Sex(ismember(tbl.ID,sub_id_exp1));
gender_exp2 = tbl.Sex(ismember(tbl.ID,sub_id_exp2));

idx = ismember(tbl.ID,sub_id_exp2);
row = find(idx==1);
tbl.ID(row(18))

%%
clc
mean(age_exp1,'omitmissing')
std(age_exp1,0,1,'omitmissing')
min(age_exp1)
max(age_exp1)
sum(gender_exp1==0)
sum(gender_exp1==1)

clc
mean(age_exp2,'omitmissing')
std(age_exp2,0,1,'omitmissing')
min(age_exp2)
max(age_exp2)
sum(gender_exp2==0)
sum(gender_exp2==1)



