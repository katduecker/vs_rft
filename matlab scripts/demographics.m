%% VS RFT
% PhD project 2

% demographics of usable participants
clear all; close all; clc
pth = 'Z:\Visual Search RFT';

addpath(fullfile(pth,'matlab scripts','alpha'))

dtpth = fullfile(pth,'data');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));
folds = subj;
clear subj

subj_age = zeros(1,length(folds));
subj_sex = cell(1,length(folds));
subj_handed_augm = zeros(1,length(folds));
for s=1:length(folds)
    matfile = dir(fullfile(dtpth,folds{s}));
    matfile = {matfile.name};
    idx = (cell2mat(cellfun(@(x) ~isempty(x), regexp(matfile,folds{s}(end-3:end)),'UniformOutput',false)) + ...
        cell2mat(cellfun(@(x) ~isempty(x), regexp(matfile,'.mat'),'UniformOutput',false))) == 2;
    
    load(fullfile(dtpth,folds{s},matfile{idx}))
    subj_age(s) = subj.age;
    subj_sex{s} = subj.gender;
    subj_handed_augm(s) = subj.augmhdsc;
end

m_age = mean(subj_age);
std_age = std(subj_age);
num_f = sum(strcmp(subj_sex,'female'));

m_hdd = mean(subj_handed_augm);
std_hdd = std(subj_handed_augm);