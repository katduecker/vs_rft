%% Pilot analysis: VS, alpha and RFT

%% Settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
ft_defaults;
subjcode = 'b59c';
soi = {'MEG2032','MEG2033','MEG2112','MEG2113','MEG2042','MEG2043','MEG1922','MEG1923','MEG2342', 'MEG2343'};
soicmb = {'MEG2032+2033','MEG2112+2113','MEG2042+2043','MEG1922+1923','MEG2342+2343'};
pthout = fullfile(mpth,'pilot','results',subjcode);

figpth = fullfile(mpth,'pilot','results','plots',subjcode);
mkdir(figpth)
fs = 1000;
load(fullfile(rsppth,[subjcode,'.mat']))
% load trigger
load(fullfile(cdpth,'trigdef.mat'))

d = dir(fullfile(megpth,subjcode));
files = {d.name};
files = files(strncmp(files,subjcode,4));

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
    trl = kd_trlfun(cfg);
    cfg.detrend = 'yes';
    cfg.trl = trl;
    %cfg.demean = 'yes';
    %cfg.baselinewindow = [-1.25 -.1];
    %cfg.hpfilter = 'yes';
    %cfg.hpfreq   = 40;
    block{fl} = ft_preprocessing(cfg);
    
end

data = ft_appenddata([],block{:})

%% Sanity checks

% plot trigger channel
trgchn = find(strcmp(data.label,'STI101'));

fig = figure;
for trl = 1:size(data.trial,2)
    curtrig = data.trial{trl}(trgchn,:);
    curtrig(curtrig > 36) = 50;
    curtrig(curtrig < 0) = 0;
    plot(data.time{trl},curtrig)
    hold on
end

% plot diode
% plot diode
fig = figure;
misc = find(strcmp(data.label,'MISC004'));
for trl = 1:size(data.trial,2)
    plot(data.time{trl},data.trial{trl}(misc,:))
    xlim([0 0.5])
    hold on
end


% artefact rejection
cfg = [];
cfg.channel = 'MEGGRAD';
datagrad = ft_selectdata(cfg,data);

cfg = [];
cfg.method = 'summary';
cfg.layout = 'neuromag306planar.lay';
datacl = ft_rejectvisual(cfg,datagrad);

[~,rej] = intersect(datagrad.sampleinfo(:,1),datacl.cfg.artfctdef.summary.artifact(:,1));

cfg = [];
keeptrl = [1:length(data.trial)];
keeptrl(rej) = [];
cfg.trials = keeptrl;
datacl = ft_selectdata(cfg,data);

% cut out long trials only
l = cell2mat(cellfun(@length,datacl.time,'UniformOutput',false));

l_trl = l >= 4500;

cfg = [];
cfg.trials = l_trl;
datal = ft_selectdata(cfg,datacl);
cfg = [];
cfg.latency = [-2.5 2];
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
datal = ft_selectdata(cfg,datal);
%% ERF

cfg = [];
cfg.keeptrials = 'no';
ERF = ft_timelockanalysis(cfg,datal);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [40 80];
%cfg.demean = 'yes';
%cfg.baselinewindow = [-1.5 0];
ERFhp = ft_preprocessing(cfg,ERF);

cfg = [];
cfg.layout = 'neuromag306planar.lay';
cfg.baseline = [-1.5 -0.1];
cfg.xlim = [-1 1];
%ft_multiplotER(cfg,ERFhp);

for c = 1:length(soi)
    plot(ERFhp.time,ERFhp.avg(strcmp(ERFhp.label,soi{c}),:))
    xlim([-1,1])
    title(soi{c})
    xlabel('time (s)')
    ylabel('amplitude fT')
    print(fullfile(figpth,[subjcode,'_',soi{c},'_ERF_bp']),'-dpng')
    pause
    close all
end


% TFR
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
cfg.foi = 0:1:100;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi));
cfg.toi = -1:0.05:.5;
cfg.keeptrials = 'no';
TFR = ft_freqanalysis(cfg,ERF);
TFR.grad = mGrad;
% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);

% plot
cfg = [];
cfg.colorbar = 'yes';
cfg.baseline = [-1.5 -0.1];
cfg.baselinetype = 'relchange';
%cfg.zlim = [-2.5 2.5];
cfg.ylim = [55 70];
cfg.xlim = [-1 0.5];
cfg.layout = 'neuromag306cmb.lay';
%ft_multiplotTFR(cfg,TFR)
% 

for c = 1:length(soicmb)
    cfg.channel = soicmb{c};
    %cfg.zlim = [0 4]
    ft_singleplotTFR(cfg,TFR)
     print(fullfile(figpth,[subjcode,'TFR_pow_',soicmb{c}]),'-dpng')
     pause
   close all
end

% sanity check: spectrum of photodiode
% check length of trials and truncate
l = cell2mat(cellfun(@length,data.time,'UniformOutput',false));
l_trl = l >= 2500;

cfg = [];
cfg.channel = 'MISC004';
cfg.trials = l_trl;
misc4 = ft_selectdata(cfg,data);
cfg = [];
cfg.latency = [-1.5 1];
misc4 = ft_selectdata(cfg,misc4);

cfg = [];
cfg.keeptrials = 'no';
misc4erf = ft_timelockanalysis(cfg,misc4);

plot(misc4erf.time,misc4erf.avg)
title('MISC004')
xlabel('time (s)')
ylabel('amplitude')
xlim([-.2 1])
print(fullfile(figpth,[subjcode,'_misc4_erf']),'-dpng')
close all

% TFR
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
%cfg.channel = {'MEGGRAD','MISC004','MISC005'};
cfg.foi = 40:1:100;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi));
cfg.toi = -1:0.05:.5;
cfg.keeptrials = 'no';
TFR4 = ft_freqanalysis(cfg,misc4erf);

% plot
imagesc(TFR4.time,TFR4.freq,squeeze(TFR4.powspctrm));
axis xy
ylabel('frequency (Hz)')
xlabel('time (s)')
title('MISC004')
ylim([40 100])
print(fullfile(figpth,[subjcode,'_misc4_tfr']),'-dpng')


%% Coherence
% generate MISC signal (as done in experiment)

scrhz = 1440;           % propixx refresh
lsig  = 2;              % maximum length of signal;
ampl = .5;
meanlum = .5;
t = linspace(0, lsig-1/fs, scrhz*lsig);
% 60 Hz signal
sig = ampl*sin(2 * pi * 60 * t)+meanlum;
%plot(t,sig)
%xlim([0 0.5])
% resample
sigrs = resample(sig,fs,scrhz);
trs = linspace(0,lsig,fs*lsig);
% hold on
% plot(trs,sigrs)
% add 2.5 baseline
misc4 = [zeros(1,2.5*fs),sigrs];
% 67 Hz signal
sig = ampl*sin(2 * pi * 67 * t)+meanlum;
sigrs = resample(sig,fs,scrhz);
misc5 = [zeros(1,2.5*fs),sigrs];

idxmisc4 = find(strcmp(datal.label,'MISC004'));
idxmisc5 = find(strcmp(datal.label,'MISC005'));


% use surrogate signals as misc signals
for t = 1:length(datal.trial)
    datal.trial{t}(idxmisc4,:) = misc4;
    datal.trial{t}(idxmisc5,:) = misc5;
end
%% hilbert transform
% to extract power spectral density of MEG sensors and diode and
% cross-spectral density between them
foi = [40:2:80];                % center frequencies of bp filter
frqwdth = 3,
cfg = [];
cfg.bpfilter    = 'yes';
cfg.hilbert     = 'complex';
cfg.keeptrials  = 'yes';

% loop over frequencies and do hilbert transform of data
for f = 1:length(foi)
    % bp filter
    cfg.bpfreq = [foi(f)-frqwdth foi(f)+frqwdth];

    data_foi = ft_preprocessing(cfg,datal);
    

    % for each channel
    for c = 1:length(datal.label)
        meg_mag = [];
        misc4_mag = [];
        
        % loop over trials
        for t = 1:length(datal.trial)
            % magnitude per trial
            meg_mag(:,t) = data_foi.trial{t}(c,:);
            misc4_mag(:,t) = data_foi.trial{t}(idxmisc4,:);
            misc5_mag(:,t) = data_foi.trial{t}(idxmisc5,:);
        end
        
        psd_meg(c,f,:) = mean(meg_mag.*conj(meg_mag),2);
        % average psd over trials -> copy such that size of psd of misc is the same
        % as psd of meg
        psd_misc4(c,f,:) = mean(misc4_mag.*conj(misc4_mag),2);
        psd_misc5(c,f,:) = mean(misc5_mag.*conj(misc5_mag),2);

        
        %csd by hand
        csd_meg_misc4(c,f,:) = abs(mean(meg_mag .* conj(misc4_mag),2)).^2;
        csd_meg_misc5(c,f,:) = abs(mean(meg_mag .* conj(misc5_mag),2)).^2;

    end
    
    
end
coh_meg_misc4_spct = csd_meg_misc4./(psd_meg.*psd_misc4);
coh_meg_misc5_spct = csd_meg_misc5./(psd_meg.*psd_misc5);
save(fullfile(pthout,[subjcode,'coh_hilbert_ftrip_fwidth_',num2str(frqwdth),'.mat']),'coh_meg_misc4_spct','coh_meg_misc5_spct','csd_meg_misc4','csd_meg_misc5','psd_meg','psd_misc4','psd_misc5')

