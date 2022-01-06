%% VS + RFT
% PhD project 2

% artefact suppression using ICA
% [c] Katharina Duecker

function a_ICA(s)

%clear all; close all; clc; beep off
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
%addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

if exist(fullfile(pth,'results','meg','3 ICA',[subjfolds{s},'_ica.mat']))
    delete(fullfile(pth,'results','meg','3 ICA',[subjfolds{s},'_ica.mat']))
end

% load in trial structure
load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'))

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-2.5*fs,meginfo.alltrl_bl{p}(:,3)+fs*2,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];
    trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;
end

% % read in artefactual sensors
% fid = fopen(fullfile(dtpth,subjfolds{s},'artef_sens.txt'),'r');
% noisy_sens = {};
% sc = 1;
% sens = fgetl(fid);
% while isstr(sens)
%     noisy_sens{sc} = sens; 
%     sc = sc+1;
%     sens = fgetl(fid);
% end
% fclose(fid)
% load in maxfiltered data using trial structure

% list fif files
d = dir(fullfile(dtpth,subjfolds{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

disp(['loading data subj ', subjfolds{s}])
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = 'MEG';
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);   
end

data = ft_appenddata([],dtprt{:});

% keep trials
keep_trl = ones(1,size(meginfo.meg_rt,1));
keep_trl(meginfo.rejtrl_all) = 0;

% discard strange trials
cfg = [];
cfg.trials = logical(keep_trl);
data = ft_selectdata(cfg,data);

% resample to 250 Hz (1/4 of sampling rate)
disp(['resample subj ', subjfolds{s}])
cfg = [];
cfg.resamplefs = 250;
dataresamp = ft_resampledata(cfg,data);

%% ICA per subject
disp(['ICA ', subjfolds{s}])
cfg = [];
cfg.method = 'runica';
cfg.numcomponent = 68;                      % identified based on first 20 subjects
dataICA = ft_componentanalysis(cfg,dataresamp);

mkdir(fullfile(pth,'results','meg','3 ICA'))
% if exist(fullfile(pth,'results','meg','3 ICA',[subjfolds{s},'_ica.mat']))
%     delete(exist(fullfile(pth,'results','meg','3 ICA',[subjfolds{s},'_ica.mat'])))
% end
save(fullfile(pth,'results','meg','3 ICA',[subjfolds{s},'_ica.mat']),'dataICA','-v7.3')
clear data* dtprt
