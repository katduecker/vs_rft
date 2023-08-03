%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a2. store matfile output from a1 into .csv files and make scatter plots
% (Supplementary Fig. 9) 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Eye movement analysis
% a: ocular artefacts for alpha high low
% b: ocular artefacts for fast vs slow trials

clear all; close all; clc; beep off;
pth = 'Z:\Visual Search RFT';

inpth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid');
mergepth = fullfile(pth,'results','meg','2 merged edf mat');
dtpth = fullfile(pth,'data');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions','alpha RFT');
occupth = fullfile(pth,'results','eyelink');
mkdir(occupth)
addpath(genpath(fullfile(pth,'matlab scripts')))
condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));

time_oi = [0.25 0.5];       % window of the median split
time_oi_str = arrayfun(@num2str, time_oi.*1000,'UniformOutput',false);


blink = zeros(length(subj),length(condi),2,2);
sacc = zeros(length(subj),length(condi),2,2);
t_bias = zeros(length(subj),length(condi),2);

for s = 1:length(subj)
    load(fullfile(occupth,subj{s},['alpha_high_low_',strjoin(time_oi_str,'_'),'.mat']))

    blink(s,:,:,:) = bl_subj_condi_bslsearch_hl;
    sacc(s,:,:,:) = sac_subj_condi_bslsearch_hl;
    t_bias(s,:,:) = bias_subj_condi_hl;

    clear bl_subj_condi_bslsearch_hl sac_subj_condi_bslsearch_hl bias_subj_condi_hl
end

%% Convert high-dim matrix into table


% bias
biasT = [];

for c = 1:length(condi)
    biasT = [biasT;[1:length(subj)]',ones(length(subj),1).*c,ones(length(subj),1),t_bias(:,c,1);[1:length(subj)]',ones(length(subj),1).*c,zeros(length(subj),1),t_bias(:,c,2)];
end


biasT = array2table(biasT,'VariableNames',{'id','condi','high','bias'});

% saccades
sacT = [];

for c = 1:length(condi)
    sacT = [sacT;[1:length(subj)]',ones(length(subj),1).*c,ones(length(subj),1),sacc(:,c,2,1);...
        [1:length(subj)]',ones(length(subj),1).*c,zeros(length(subj),1),sacc(:,c,2,2)];
end

sacT = array2table(sacT,'VariableNames',{'id','condi','high','num_sac'});


% blinks
blinkT = [];

for c = 1:length(condi)
    blinkT = [blinkT;[1:length(subj)]',ones(length(subj),1).*c,ones(length(subj),1),blink(:,c,2,1);...
        [1:length(subj)]',ones(length(subj),1).*c,zeros(length(subj),1),blink(:,c,2,2)];
end

blinkT = array2table(blinkT,'VariableNames',{'id','condi','high','num_bl'});



%% Plots

fig = figure('Position',[0 0 1920 1080/2.5]);
subplot(141)
scatter(t_bias(:,[2,1,4,3],1),t_bias(:,[2,1,4,3],2),[],col_palette,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
scatter(mean(t_bias(:,[2,1,4,3],1),1),mean(t_bias(:,[2,1,4,3],2),1),50,col_palette,'filled','MarkerEdgeColor',[0 0 0])
diag = linspace(.35,.65,35);
xlim([0.35 0.65])
ylim([0.35 0.65])
xticks(0.35:0.15:0.65)
yticks(0.35:0.15:0.65)
xtickformat('%.2f')
ytickformat('%.2f')
plot(diag,diag,'-.','Color','k')
pbaspect([1 1 1])
xlabel('target bias fast')
ylabel('target bias slow')
print(fig,fullfile(occupth,['target_bias_alpha_high_low',strjoin(time_oi_str,'_')]),'-dsvg')
print(fig,fullfile(occupth,['target_alpha_high_low',strjoin(time_oi_str,'_')]),'-dpng')

% saccades and blinks
fig = figure('Position',[0 0 1920 1080/2.5]);


subplot(141)
x = [blink(:,[2,1,4,3],2,1);blink(:,[2,1,4,3],2,2)];
mx = round(max(x(:))+0.05,1);
scatter(blink(:,[2,1,4,3],2,1),blink(:,[2,1,4,3],2,2),[],col_palette,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
scatter(mean(blink(:,[2,1,4,3],2,1),1),mean(blink(:,[2,1,4,3],2,2),1),50,col_palette,'filled','MarkerEdgeColor',[0 0 0])
diag = linspace(0,mx,140);
xticks(0:mx/2:mx)
yticks(0:mx/2:mx)
xtickformat('%.1f')
ytickformat('%.1f')
plot(diag,diag,'-.','Color','k')
pbaspect([1 1 1])
xlabel('blinks/trial fast')
ylabel('blinks/trial slow')

subplot(142)
x = [sacc(:,[2,1,4,3],2,1);sacc(:,[2,1,4,3],2,2)];
mx = round(max(x(:)));
scatter(sacc(:,[2,1,4,3],2,1),sacc(:,[2,1,4,3],2,2),[],col_palette,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
hold on
scatter(mean(sacc(:,[2,1,4,3],2,1),1),mean(sacc(:,[2,1,4,3],2,2),1),50,col_palette,'filled','MarkerEdgeColor',[0 0 0])
diag = linspace(0,mx,140);
xticks(0:mx/2:mx)
yticks(0:mx/2:mx)
xlim([0 mx])
ylim([0 mx])
xtickformat('%.1f')
ytickformat('%.1f')
plot(diag,diag,'-.','Color','k')
pbaspect([1 1 1])
xlabel('saccades/trial fast')
ylabel('saccades/trial slow')
%title('search')

print(fig,fullfile(occupth,['artefacts_alpha_high_low_',strjoin(time_oi_str,'_')]),'-dsvg')
print(fig,fullfile(occupth,['artefacts_alpha_high_low_',strjoin(time_oi_str,'_')]),'-dpng')



writetable(biasT,fullfile(occupth,['bias_alpha_high_low_',strjoin(time_oi_str,'_'),'.csv']))
writetable(sacT,fullfile(occupth,['saccades_alpha_high_low_',strjoin(time_oi_str,'_'),'.csv']))
writetable(blinkT,fullfile(occupth,['blinks_alpha_high_low_',strjoin(time_oi_str,'_'),'.csv']))
