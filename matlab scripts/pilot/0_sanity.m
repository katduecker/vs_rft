%% Sanity check scripts Visual Search

clear all; close all; clc; beep off;

% settings
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'data');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts');
addpath(genpath(mtlpth))
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
%addpath('Z:\fieldtrip')
ft_defaults;
subjcode = 'b3ec';

d = dir(megpth)
folds = {d.name};
folds = folds(strncmp(folds,'202',2))
f = cell2mat(cellfun(@(x) strcmp(x(end-3:end),subjcode),folds,'UniformOutput',false));
subpth = fullfile(megpth,folds{f},'meg');
d = dir(subpth)
files = {d.name};
files = files(strncmp(files,'part',4));

clear d f folds

% average grad positions

% load in grad structures
grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(subpth,files{fl}))];
end
% average
mGrad = grad(1);
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
clear grad

%% Load in data

% load trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(subpth,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    % target or distractor tagging?
    cfg.detrend = 'yes';
    
    % load in all trials
    cfg.trl = kd_trlfun_alltrl(cfg); 
    block{fl} = ft_preprocessing(cfg);
end

data = ft_appenddata([],block{:});