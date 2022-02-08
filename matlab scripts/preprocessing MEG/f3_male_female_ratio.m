%% VS + RFT
% PhD project 2

% Male/female ratio
clear all; close all; clc; beep off;
pth = 'Z:\Visual Search RFT';
dtpth = fullfile(pth,'data');
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');

% list usable subjects
d = dir(dtpth);
d = {d.name};
subjfolds = d(strncmp(d,'202',3));

load(fullfile(scriptpth,'idx_subjoi.mat'))

subjfolds = subjfolds(usable_idx);
clear d

subj_sx = [];

for s = 1:length(subjfolds)
    % load response file

    d = dir(fullfile(dtpth,subjfolds{s}));
    d = {d.name};
    
    % find indices with subj id
    id = cell2mat(cellfun(@(x) ~isempty(x), regexp(d,subjfolds{s}(end-3:end)),'UniformOutput',false));
    d = d(id);
    id = cell2mat(cellfun(@(x) ~isempty(x),regexp(d,'mat'),'UniformOutput',false));
    
    load(fullfile(dtpth,subjfolds{s},d{id}))

    if  strcmp(subj.gender,'female')
        subj_sx = [subj_sx, 1];
    elseif  strcmp(subj.gender,'male')
        subj_sx = [subj_sx,2];
    else
        subj_sx = [subj_sx,0];
    end

    clear subj
    
end

numf = length(find(subj_sx == 1))
numm = length(find(subj_sx == 2))


