%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b2: plot RT and hit-rate difference for alpha high vs low (Fig 4. c & h)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Behavioural analyses
% a: performance per condition
% b. performance for alpha high vs low


%% Paths
clear all; close all; clc

split_ta_tp = 1;            % split for target present/absent?
toi_alpha_split = [0.25 0.5];   % alpha toi
bslcor = 1;   
set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Arial')


addpath('Z:\fieldtrip')
ft_defaults;

pth = 'Z:\Visual Search RFT';
addpath('Z:\Visual Search RFT\Violinplot-Matlab-master')
addpath(genpath('Z:\Visual Search RFT\ScientificColourMaps7\'))
load('berlin.mat')
behavpth = fullfile(pth,'results','behavior');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','alpha high low');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;
toi_alpha_split = arrayfun(@num2str,toi_alpha_split.*1000,'UniformOutput',false);


if split_ta_tp

    condi = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},{'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};
    split_suf = '_ta_tp';

else
    condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
    split_suf = '';

end

varnames = {'ung 16','gui 16', 'ung 32','gui 32'};

load(fullfile(behavpth,['RT_hitrate_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.mat']))

fig = figure('Position',[0 0 1920/4 1080/1.5]);
subplot(121)
ax = violinplot(avg_rt_diff,{'RT'},'ViolinColor',[0 0 0],'ShowMean',true);
ylim([-0.3 0.3])
yticks([-0.3:0.3:0.3])

subplot(122)
ax = violinplot(avg_hit_diff,{'hit'},'ViolinColor',berlin(end-50,:),'ShowMean',true);
ylim([-0.1 0.1])
yticks([-0.1:0.1:0.1])

print(fig,fullfile(behavpth,['RT_hitrate_diff_avg',strjoin(toi_alpha_split,'_'),split_suf]),'-dpng','-r0')
print(fig,fullfile(behavpth,['RT_hitrate_diff_avg',strjoin(toi_alpha_split,'_'),split_suf]),'-dsvg','-r0')
