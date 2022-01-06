%% VS + RFT
% PhD project 2

% check which participants have to be excluded based on the soi
% identification
clear all; close all; clc
pth = 'Z:\Visual Search RFT';
addpath('Z:\fieldtrip')
% addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb','soi');

addpath(fullfile('W:\','fieldtrip'))
ft_defaults;

% list subj
d = dir(cohpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

load(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))

empty_idx = logical(cell2mat(cellfun(@isempty,mergesubj(:,2),'UniformOutput',false)));
mergesubj = mergesubj(~empty_idx,:)
mergesubj = [mergesubj,cell(length(subjfolds),1)];

for s = 1:length(subjfolds)
    try 
        load(fullfile(cohpth,subjfolds{s},'soi_stat.mat'))
    if isempty(soi_stat)
        mergesubj{s,5} = 'no soi';
    else
        mergesubj{s,5} = 'soi identified';
    end
    clear soi_stat
    catch ME
    end
end
save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_soi_coh.mat'),'mergesubj')

numberusable = sum(strcmp(mergesubj(:,5),'soi identified'));

% list all p's with raw data
d = dir(megpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
[~,usable_idx] = intersect(subjfolds,mergesubj(strcmp(mergesubj(:,5),'soi identified'),1))

save(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'),'usable_idx')