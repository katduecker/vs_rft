%% VS + RFT
% PhD project 2

% semi-automatic artefact rejection

% [c] Katharina Duecker

function a2_fun_artef_rej(s,pth, fs)

dtpth = fullfile(pth,'data');                                       % raw data
maxfpth = fullfile(pth,'results','meg', '1 maxfilter');             % max filter
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');   % scripts
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
load(fullfile(scriptpth,'idx_subjoi.mat'))                          % usable subjects index

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

subj = subjfolds(usable_idx);

% load in trial structure
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*1.5,meginfo.alltrl_bl{p}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];
    trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;
end

% list maxfiltered data
d = dir(fullfile(maxfpth,subj{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load data
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(maxfpth,subj{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = 'MEG';
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);   
end

data = ft_appenddata([],dtprt{:});

% discard strange trials
cfg = [];
cfg.trials = ~ismember(1:length(data.trial),meginfo.rejtrl_all);
data = ft_selectdata(cfg,data);

% sampleinfo
trl_end_samp = cellfun(@length,data.trial,'UniformOutput',false).*[1:length(data.trial)];
data.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];

%% Artefact identification

% split into grad and mag
cfg = [];
cfg.channel = 'MEGGRAD';
graddata = ft_selectdata(cfg,data);
cfg.channel = 'MEGMAG';
magdata = ft_selectdata(cfg,data);


% reject visual
cfg = [];
cfg.method = 'summary';
cfg.layout = 'neuromag306planar.lay';
grad_rej = ft_rejectvisual(cfg,graddata);

% get idx of to be kept trials
[~, keep_trl_grad] = intersect(graddata.sampleinfo(:,1),grad_rej.sampleinfo(:,1));

cfg.layout = 'neuromag306mag.lay';
mag_rej = ft_rejectvisual(cfg,magdata);
[~, keep_trl_mag] = intersect(magdata.sampleinfo(:,1),mag_rej.sampleinfo(:,1));

keep_trl_idx = intersect(keep_trl_grad,keep_trl_mag);

keep_trl = zeros(1,length(data.trial));
keep_trl(keep_trl_idx) = 1;
find(~keep_trl)

meginfo.post_artefrej_keep = keep_trl;
save(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'),'meginfo','-append')