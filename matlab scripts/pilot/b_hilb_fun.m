%% Pilot analysis: VS, alpha and RFT
% coherence using hilbert transform as function

function b_hilb_fun(frqwdth)

%% Settings
mpth = '/rds/projects/j/jenseno-visual-search-rft';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
ft_defaults;
subjcode = 'b57a';


soi = {'MEG2032','MEG2033','MEG2112','MEG2113'};
soicmb = {'MEG2032+2033','MEG2112+2113','MEG2042+2043','MEG1922+1923','MEG2342+2343'};
pthout = fullfile(mpth,'pilot','results',subjcode);

figpth = fullfile(mpth,'pilot','results','plots',subjcode);
mkdir(figpth)
fs = 1000;
load(fullfile(rsppth,[subjcode,'.mat']))

load(fullfile(cdpth,'trigdef.mat'))

d = dir(fullfile(megpth,subjcode));
files = {d.name};
files = files(strncmp(files,subjcode,4));

% just check color for now
files(cell2mat(cellfun(@isempty,strfind(files,'col'),'UniformOutput',false))) = [];

%% Load in data
% separately: misc004 picks up 60 Hz signal and misc005 picks up 60 Hz
% signal
%event = [];
grad = [];
for fl = 1:length(files)
   % event = [event;ft_read_event(fullfile(megpth,subjcode,files{fl}))];
    grad = [grad;ft_read_sens(fullfile(megpth,subjcode,files{fl}))];
end

% mean channel position
mGrad = grad(1);
% average grad structure
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);


% read in trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    [trlsmp_misc004{fl}, trltrg_misc004{fl}, trl_misc004{fl},trlsmp_misc005{fl}, trltrg_misc005{fl}, trl_misc005{fl}] ...
        = kd_trlfun_trl_phd(cfg);
    cfg.trl = trl_misc004{fl};
    cfg.detrend = 'yes';
    %cfg.demean = 'yes';
    %cfg.baselinewindow = [-1.25 -.1];
    %cfg.hpfilter = 'yes';
    %cfg.hpfreq   = 40;
    block_misc004{fl} = ft_preprocessing(cfg);
    cfg.trl = trl_misc005{fl};
    block_misc005{fl} = ft_preprocessing(cfg);

end

% 1 block empty?
block_misc004(cell2mat(cellfun(@(x) isempty(x.trial),block_misc004,'UniformOutput',false)))= [];
block_misc005(cell2mat(cellfun(@(x) isempty(x.trial),block_misc005,'UniformOutput',false)))= [];

% concatenate
cfg = [];
cfg.keepsampleinfo = 'no';
data.misc004 = ft_appenddata(cfg,block_misc004{:});
data.misc005 = ft_appenddata(cfg,block_misc005{:}); 

% cut out long trials
l_004 = cell2mat(cellfun(@length,data.misc004.time,'UniformOutput',false));
l_005 = cell2mat(cellfun(@length,data.misc005.time,'UniformOutput',false));

% select data trials >= 1s (+2 second before and after
cfg = [];
cfg.trials = l_004 >= 4500;
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
data.misc004l = ft_selectdata(cfg,data.misc004);
cfg.trials = l_005 >= 4500;
data.misc005l = ft_selectdata(cfg,data.misc005);

cfg = [];
cfg.latency = [-2.5 2];
data.misc004l = ft_selectdata(cfg,data.misc004l);
data.misc005l = ft_selectdata(cfg,data.misc005l);

%% hilbert transform
% to extract power spectral density of MEG sensors and diode and
% cross-spectral density between them
foi = [40:2:80];                % center frequencies of bp filter

cfg = [];
cfg.bpfilter    = 'yes';
cfg.hilbert     = 'complex';
cfg.keeptrials  = 'yes';
idxmisc004 = find(strcmp(data.misc004l.label,'MISC004'));

% loop over frequencies and do hilbert transform of data
for f = 1:length(foi)
    % bp filter
    cfg.bpfreq = [foi(f)-frqwdth foi(f)+frqwdth];

    data_foi = ft_preprocessing(cfg,data.misc004l);
    

    % for each channel
    for c = 1:length(data.misc004l.label)
        meg_mag = [];
        misc4_mag = [];
        
        % loop over trials
        for t = 1:length(data.misc004l.trial)
            % magnitude per trial
            meg_mag(:,t) = data_foi.trial{t}(c,:);
            misc4_mag(:,t) = data_foi.trial{t}(idxmisc004,:);

        end
        psd_meg(c,f,:) = mean(meg_mag.*conj(meg_mag),2);
        % average psd over trials -> copy such that size of psd of misc is the same
        % as psd of meg
        psd_misc4(c,f,:) = mean(misc4_mag.*conj(misc4_mag),2);
        
        
        %csd by hand
        csd_meg_misc4(c,f,:) = abs(mean(meg_mag .* conj(misc4_mag),2)).^2;
        
    end
    
    
end

coh_meg_misc4_spct = csd_meg_misc4./(psd_meg.*psd_misc4);

% clean up workspace
clear data* l_* *_mag trl*

% load template coherence
load(fullfile(mpth,'pilot','results',[subjcode,'_COH_misc004.mat']))

coh_meg_misc4 = COH004;
clear COH004

% add respective fields to template structure
coh_meg_misc4.cohspctrm = coh_meg_misc4_spct;
coh_meg_misc4.freq = foi;
coh_meg_misc4.time = data.misc004l.time{1};
coh_meg_misc4.dof  = repmat(t,1,length(foi));
coh_meg_misc4.grad = mGrad;
coh_meg_misc4 = rmfield(coh_meg_misc4,'cfg');
coh_meg_misc4.labelcmb = [coh_meg_misc4.labelcmb;{'MISC004'}, {'MISC004'};{'MISC005'}, {'MISC005'}];
coh_meg_misc4.label = coh_meg_misc4.labelcmb(:,1);
save(fullfile(pthout,[subjcode,'coh_hilbert_ftrip_fwidth_',num2str(frqwdth),'.mat']),'coh_meg_misc4','csd_meg_misc4','psd_meg','psd_misc4')
