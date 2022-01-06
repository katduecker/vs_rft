%% Pilot analysis: VS, alpha and RFT

% coherence using hilbert transform
% Coherence part (c) Y.Pan, K.Duecker

%% Settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
pthout = fullfile(mpth,'pilot','results');
addpath(mtlpth)
ft_defaults;
subjcode = 'b57a';
figpth = fullfile(mpth,'pilot','results','plots',subjcode);
mkdir(figpth)
fs = 1000;
load(fullfile(rsppth,[subjcode,'.mat']))
% load trigger
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

% separate responses accordingly
[rspns.misc004, rspns.misc005] = kd_rspns_phd(subj);

% cut out long trials
l_004 = cell2mat(cellfun(@length,data.misc004.time,'UniformOutput',false));
l_005 = cell2mat(cellfun(@length,data.misc005.time,'UniformOutput',false));

% select data trials >= 1s (+2 second before and after
cfg = [];
cfg.trials = l_004 >= 4500;
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
%cfg.latency = [-1.5 1.5];
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
    cfg.bpfreq = [foi(f)-5 foi(f)+5];

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
load(fullfile(mpth,'pilot','results',[subjcode,'_COH_misc004.mat']),'COH004')

coh_meg_misc4 = COH004;
clear COH004

% add respective fields to template structure
coh_meg_misc4.cohspctrm = coh_meg_misc4_spct;
coh_meg_misc4.freq = foi;
coh_meg_misc4.time = linspace(-1.5,1,size(coh_meg_misc4_spct,3));
coh_meg_misc4.dof  = repmat(t,1,length(foi));
coh_meg_misc4.grad = mGrad;
coh_meg_misc4 = rmfield(coh_meg_misc4,'cfg');
coh_meg_misc4.labelcmb = [coh_meg_misc4.labelcmb;{'MISC004'}, {'MISC004'};{'MISC005'}, {'MISC005'}];
coh_meg_misc4.label = coh_meg_misc4.labelcmb(:,1);
save(fullfile(pthout,[subjcode,'coh_hilbert_ftrip.mat']),'coh_meg_misc4','csd_meg_misc4','psd_meg','psd_misc4')

% plot
cfg = [];
cfg.parameter = 'cohspctrm';
%cfg.channel = 'MEGGRAD';
cfg.refchannel = 'MISC004';
cfg.layout = 'neuromag306planar.lay';
cfg.colorbar = 'yes';
ft_multiplotTFR(cfg,coh_meg_misc4)
ft_hastoolbox('brewermap',1);
colormap(flipud(brewermap(64,'RdBu')))

imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(coh_meg_misc4.cohspctrm(end,:,:)))
axis xy
  

close all
plot(coh_meg_misc4.time,squeeze(psd_meg(1,end,:)))
    
    
%% Hard code without fieldtrip

foi = [52:1:75];                % center frequencies of bp filter

idxmisc004 = find(strcmp(data.misc004l.label,'MISC004'));
idxmisc005 = find(strcmp(data.misc005l.label,'MISC005'));

% concatenate data for filtering
catdata = horzcat(data.misc004l.trial{:});

% length of one trial
trl_l = length(data.misc004l.trial{1});
% nyquist frequency:
nyqf = 1000/2;
for c = 1:length(data.misc004l.label)
    % loop over frequencies and do hilbert transform of data
    for f = 1:length(foi)
        % bp filter: butterworth
        [bbp, abp] = butter(2,[(foi(f)-2)/nyqf (foi(f)+2)/nyqf]);
        
        fdat = filter(bbp,abp,catdata(c,:));
        % this is the same thing every iteration - just to match size of
        % psd matrices
        fdatmisc4 = filter(bbp,abp,catdata(idxmisc004,:));
        fdatmisc5 = filter(bbp,abp,catdata(idxmisc005,:));
        
        fdatmisc = 
        % separate back into single trials and compute magnitude
        t = 1;          % trial counter
        ts = 1;
        while ts < size(fdat,2)
            % hilbert transform per trial
            hilbtrl(t,:) = hilbert(fdat(ts:ts+trl_l-1));
            hilbtrlmisc4(t,:) = hilbert(fdatmisc4(ts:ts+trl_l-1));
            hilbtrlmisc5(t,:) = hilbert(fdatmisc5(ts:ts+trl_l-1));
            
            
            ts = ts + trl_l;
            t = t + 1;
        end
        
        % compute psd
        psd_meg(c,f,:) = mean(hilbtrl.*conj(hilbtrl),1);
        psd_misc4(c,f,:) = mean(hilbtrlmisc4.*conj(hilbtrlmisc4),1);
        psd_misc5(c,f,:) = mean(hilbtrlmisc5.*conj(hilbtrlmisc5),1);

        % psd misc4 
        % compute |csd|²
        csd_meg_misc4(c,f,:) = abs(mean(hilbtrl.*conj(hilbtrlmisc4),1)).^2;
        csd_meg_misc5(c,f,:) = abs(mean(hilbtrl.*conj(hilbtrlmisc5),1)).^2;
        clear hilbtrl* *hilb fdat*
    end
end

% coherence misc4 60 Hz
coh_meg_misc4_spct = csd_meg_misc4./(psd_meg.*psd_misc4);
% coherence misc5 60 Hz
coh_meg_misc5_spct = csd_meg_misc5./(psd_meg.*psd_misc5);
% load template coherence
load(fullfile(mpth,'pilot','results',[subjcode,'_COH_misc004.mat']),'COH004')

coh_meg_misc4 = COH004;
clear COH004

% add respective fields to template structure
coh_meg_misc4.cohspctrm = coh_meg_misc4_spct;
coh_meg_misc4.freq = foi;
coh_meg_misc4.time = linspace(-1.5,1,size(coh_meg_misc4_spct,3));
coh_meg_misc4.dof  = repmat(t,1,length(foi));
coh_meg_misc4.grad = mGrad;
coh_meg_misc4 = rmfield(coh_meg_misc4,'cfg');
coh_meg_misc4.labelcmb = [coh_meg_misc4.labelcmb;{'MISC004'}, {'MISC004'};{'MISC005'}, {'MISC005'}];
coh_meg_misc4.label = coh_meg_misc4.labelcmb(:,1);

% check psd of misc 4
subplot(311)
imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(psd_misc4(1,:,:)))
axis xy
title('misc004: 60 Hz signal')
subplot(312)
imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(psd_misc5(1,:,:)))
axis xy
title('misc005: 67 Hz signal')
subplot(313)
% coherence misc4 misc5
imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(coh_meg_misc4.cohspctrm(idxmisc005,:,:)))
axis xy
colorbar
title('coherence misc004 misc005')

% plot
cfg = [];
cfg.parameter = 'cohspctrm';
%cfg.channel = 'MEGGRAD';
cfg.refchannel = 'MISC004';
cfg.layout = 'neuromag306planar.lay';
cfg.colorbar = 'yes';
ft_multiplotTFR(cfg,coh_meg_misc4)
ft_hastoolbox('brewermap',1);
colormap(flipud(brewermap(64,'RdBu')))

close all
soi = {'MEG2533','MEG2242','MEG0213','MEG2133','MEG2522','MEG1943'};

for s = 1:length(soi)
    cfg.channel = soi{s};
    ft_singleplotTFR(cfg,coh_meg_misc4)
    ft_hastoolbox('brewermap',1);
    colormap(flipud(brewermap(64,'RdBu')))
        pause

    print(fullfile(figpth,[subjcode,'_coh_',soi{s}]),'-dpng')
    close all
end
save(fullfile(pthout,[subjcode,'coh_hilbert_hardcode.mat']),'coh_meg_misc4','coh_meg_misc5_spct','csd_meg_misc4','csd_meg_misc5','psd_meg','psd_misc4','psd_misc5')