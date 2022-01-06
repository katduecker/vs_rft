%% VS + RFT

clear all; close all; clc
% exploration RT data of participants to be kept

pth = 'E:\UoB\Proj2_Visual Search';
rawdt = fullfile(pth, 'data');
soipth = fullfile(pth,'results','meg','5 COH hilb','soi');
mtlbdir = 'C:\Users\katha\Documents\MATLAB';
rt_pre_pth = fullfile(pth,'results','behavior','preprocessing','explore');
mkdir(rt_pre_pth)

d = dir(soipth);
fl = {d.name};
fl(1:2) = [];

%% load in soi and create logical index for subjects to be kept
keepsubj = zeros(1,length(fl));

for s = 1:length(fl)
    load(fullfile(soipth,fl{s},'soi_stat.mat'))
    
    if ~isempty(soi_stat)
        keepsubj(s) = 1;
    end
    clear soi_stat
end

subjfolds = fl(logical(keepsubj));

clear d fl keepsubj
%% create reaction time matrix for all p
load(fullfile(pth,'experiment','trigdef.mat'))
% columns: set size (2), target absent/target present (0/1), hit/miss/fa (1/2/3), RT
rt_all = nan(4,length(subjfolds),960);

for s = 1:length(subjfolds)
    load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
        'trl_overlap_meg_el_rsp.mat'),'rspinfo')
    subjcode = subjfolds{s}(end-3:end);
    d = dir(fullfile(rawdt,subjfolds{s}));
    fl = {d.name};
    % find mat files (that's converted edf and response)
    fl_idx = cell2mat(cellfun(@(x) ~isempty(x),strfind(fl,'mat'),'UniformOutput',false));
    fl = fl(fl_idx);
    % find response file
    fl_idx = cell2mat(cellfun(@(x) ~isempty(x),strfind(fl,subjcode),'UniformOutput',false));
    if length(fl(fl_idx)) > 1
        % delete practice file
        fl = fl(cell2mat(cellfun(@isempty,strfind(fl(fl_idx),'pract'),'UniformOutput',false)));
        
        if length(fl) > 1
            error('subject has more than one response file, merge before preproc')
        end
        
        load(fullfile(rawdt,subjfolds{s},fl{1}));

    else
        load(fullfile(rawdt,subjfolds{s},fl{fl_idx}));
    end

    % trial specifics, transpose for linear indexing (see merge script)
    trl_specs = [];
    for b = 1:size(subj.exp.trials{1},1)
        trl_specs = [trl_specs,subj.exp.trials{1}(b,:)];
    end
    
    %triple check if trl_specs indexing and rspinfo are the same
    comp_rsp_trl = zeros(960,1);
    for t = 1:prod(size(trl_specs))
        curtrl = cellfun(@num2str,trl_specs{t},'UniformOutput',false);
        curtrl = [curtrl{:}];
        comp_rsp_trl(t) = strcmp(curtrl,trigdef{rspinfo.trl{t,1},2});
        clear curtrl
    end
    
    if ~isempty(find(comp_rsp_trl == 0))
        error('trials are not the same')
    end
    
    for t = 1:prod(size(trl_specs))
        % set size
        rt_all(1,s,t) = trl_specs{t}{2};
        % target present/absent
        rt_all(2,s,t) = strcmp(trl_specs{t}{3},'tp');
    end
    
    % code hits
    rt_all(3,s,strcmp(rspinfo.trl(:,2),'h')) = 1;
    % code miss
    rt_all(3,s,strcmp(rspinfo.trl(:,2),'m')) = 2;
    % code false alarm
    rt_all(3,s,strcmp(rspinfo.trl(:,2),'fa')) = 3;
    
    % reaction time 
    rt_all(4,s,:) = [rspinfo.trl{:,3}];
    
    clear rspinfo trl_specs comp_rsp_trl
end

%% Raincloud plot
[cb] = cbrewer('qual', 'Set3', length(subjfolds), 'pchip');


addpath(fullfile(mtlbdir,'cbrewer'))
addpath(fullfile(mtlbdir,'RainCloudPlots','tutorial_matlab'))

fig = figure('Position',[0 0 800 900]);
% all subjects
subplot(311)
h1 = raincloud_plot(reshape(rt_all(4,:,:),[],1), 'box_on',1);
title('all subjects all conditions')

% set size
subplot(312)
idx_16 = rt_all(1,:,:) == 16;
h1 = raincloud_plot(reshape(rt_all(4,squeeze(idx_16)),[],1), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
 
h2 = raincloud_plot(reshape(rt_all(4,squeeze(~idx_16)),[],1), 'box_on',1,'color', cb(4,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
     'box_col_match', 1);
legend([h1{1} h2{1}], {'set size 16', 'set size 32'});
title('all subjects per set size')

% target present/absent
subplot(313)
tp = logical(squeeze(rt_all(2,:,:)));
h1 = raincloud_plot(reshape(rt_all(4,tp),[],1), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
 
h2 = raincloud_plot(reshape(rt_all(4,~tp),[],1), 'box_on',1,'color', cb(3,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
     'box_col_match', 1);
 
legend([h1{1} h2{1}], {'target present', 'target absent'});
title('all subjects target present vs absent')
print(fig,fullfile(rt_pre_pth,'rt_explore_condi'),'-dpng')

% per subject
fig = figure('Position',[0 0 800 900]);
for s = 1:length(subjfolds)
    h1 = raincloud_plot(reshape(rt_all(4,s,:),[],1), 'box_on', 0, 'color', cb(s,:), 'alpha', 0.5);
end
title('all subjects target present vs absent')
print(fig,fullfile(rt_pre_pth,'rt_explore_condi'),'-dpng')


fig = figure('Position',[0 0 1600 900]);

s = 0;

for b = 1:4
    subplot(2,2,b)
    for x = 1:4
        s = s +1;
        h1 = raincloud_plot(reshape(rt_all(4,s,:),[],1), 'box_on', 1, 'color', cb(s,:), 'alpha', 0.5,...
            'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
            'box_col_match', 1);
    end
end
print(fig,fullfile(rt_pre_pth,'rt_explore_subjects_split'),'-dpng')
