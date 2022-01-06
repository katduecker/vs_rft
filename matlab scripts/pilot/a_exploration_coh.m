%% Pilot analysis: VS, alpha and RFT

%% Settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
ft_defaults;
subjcode = 'b58b';

pthout = fullfile(mpth,'pilot','results',subjcode);

figpth = fullfile(mpth,'pilot','results','plots',subjcode);
mkdir(figpth)
fs = 1000;

% load trigger
load(fullfile(cdpth,'trigdef.mat'))

d = dir(fullfile(megpth,subjcode));
files = {d.name};
files = files(strncmp(files,'part',4));


%% Load in data

% average grad positions
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
clear grad

%% Separate into 60 and 67 Hz stimulation

tagblock = {'t','t','d','d'};
% load trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    
    % target or distractor tagging?
    cfg.tag = tagblock{fl};
    [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg);
%     cfg.hpfilter = 'yes';
%     cfg.hpfreq   = 52;
%     cfg.detrend = 'yes';
    
    % target tagged at 60 Hz
    cfg.trl = trl60_misc4;  
    block60misc4{fl} = ft_preprocessing(cfg);
    % target tagged at 60 Hz
    cfg.trl = trl60_misc5;  
    block60misc5{fl} = ft_preprocessing(cfg);
    
     % target tagged at 67 Hz
    cfg.trl = trl67_misc4;  
    block67misc4{fl} = ft_preprocessing(cfg);  
    
    % target tagged at 67 Hz
    cfg.trl = trl67_misc5;  
    block67misc5{fl} = ft_preprocessing(cfg);  
end

data60misc4 = ft_appenddata([],block60misc4{:});
data60misc5 = ft_appenddata([],block60misc5{:});

data67misc4 = ft_appenddata([],block67misc4{:});
data67misc5 = ft_appenddata([],block67misc5{:});


% discard strange trials
l = cell2mat(cellfun(@length,data60misc4.time,'UniformOutput',false));
cfg = [];
cfg.channel = {'MEGGRAD','MISC004'};
cfg.trials = l./fs < 10;
data60misc4 = ft_selectdata(cfg,data60misc4);
% discard strange trials
l = cell2mat(cellfun(@length,data60misc5.time,'UniformOutput',false));
cfg.trials = l./fs < 10;
cfg.channel = {'MEGGRAD','MISC005'};
data60misc5 = ft_selectdata(cfg,data60misc5);

% discard strange trials
l = cell2mat(cellfun(@length,data67misc4.time,'UniformOutput',false));
cfg.trials = l./fs < 10;
cfg.channel = {'MEGGRAD','MISC004'};

data67misc4 = ft_selectdata(cfg,data67misc4);

% discard strange trials
l = cell2mat(cellfun(@length,data67misc5.time,'UniformOutput',false));
cfg.trials = l./fs < 10;
cfg.channel = {'MEGGRAD','MISC005'};

data67misc5 = ft_selectdata(cfg,data67misc5);


%% Calculate coherence by hand

foi = 52:2:80;
frqwdth = 2;                  % width passband          
[coh60misc4_spect, psdmeg60misc4, psdmisc60misc4] = kd_coh_hilb(data60misc4, foi,frqwdth, 'MISC004');
[coh60misc5_spect, psdmeg60misc5, psdmisc60misc5] = kd_coh_hilb(data60misc5, foi,frqwdth, 'MISC005');
[coh67misc4_spect, psdmeg67misc4, psdmisc67misc4] = kd_coh_hilb(data67misc4, foi,frqwdth, 'MISC004');
[coh67misc5_spect, psdmeg67misc5, psdmisc67misc5] = kd_coh_hilb(data67misc5, foi,frqwdth, 'MISC005');

save(fullfile(mpth,'pilot','results',[subjcode,'coh_tag_onegroup.mat']),'coh60misc4_spect','coh60misc5_spect','coh67misc4_spect','coh67misc5_spect')
% load template coherence struct
load(fullfile(mpth,'pilot','results','b57a_COH_misc004.mat'))

fig = figure;
% sanity check
subplot(221)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(psdmisc60misc4(1,:,:)))
title('PSD diode 1 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
subplot(222)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(psdmisc67misc4(1,:,:)))
title('PSD diode 1 67 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
subplot(223)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(psdmisc60misc5(1,:,:)))
title('PSD diode 2 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
subplot(224)
imagesc(linspace(-2.5,5.5,8000),foi,squeeze(psdmisc67misc5(1,:,:)))
title('PSD diode 2 67 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
xlim([-1.5 4])
axis xy
colorbar
print(fig, fullfile(figpth,'psd_diode'),'-dpng')




coh60misc4 = COH004;
clear COH004

% add respective fields to template structure
coh60misc4.cohspctrm = coh60misc4_spect;
coh60misc4.freq = foi;
coh60misc4.time = linspace(-2.5,5.5,8000);
coh60misc4.dof  = repmat(160,1,length(foi));
coh60misc4.grad = mGrad;
coh60misc4 = rmfield(coh60misc4,'cfg');
coh60misc4.labelcmb = [coh60misc4.labelcmb;{'MISC004'}, {'MISC004'}];
coh60misc4.label = coh60misc4.labelcmb(:,1);

coh60misc5 = coh60misc4;
coh60misc5.cohspctrm = coh60misc5_spect;
coh67misc4 = coh60misc4;
coh67misc4.cohspctrm = coh67misc4_spect;
coh67misc5 = coh60misc4;
coh67misc5.cohspctrm = coh67misc5_spect;

cfg = [];
cfg.layout = 'neuromag306planar.lay';
cfg.parameter = 'cohspctrm';
cfg.refchannel = 'MISC004';
% cfg.ylim = [40 80];
% cfg.baseline = [-1.5 -0.5];
% cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.xlim = [-1.5 2];
fig = figure;
ft_multiplotTFR(cfg,coh60misc4)
print(fig,fullfile(figpth,'coh60misc4'),'-dpng')
close all


fig = figure;
%subplot SOI

soi = {'MEG2032', 'MEG2033', 'MEG2042', 'MEG2043;', 'MEG2113', 'MEG2112'};
[~,soiidx] = intersect(coh60misc4.label,soi)
subplot(221)
imagesc(coh60misc4.time,coh60misc4.freq,squeeze(mean(coh60misc4.cohspctrm(soiidx,:,:),1)));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('coherence diode 1 60 Hz')
xlim([-1.5 1])
colorbar
axis xy
caxis([0 0.02])
subplot(222)
imagesc(coh67misc4.time,coh67misc4.freq,squeeze(mean(coh67misc4.cohspctrm(soiidx,:,:),1)));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('coherence diode 1 67 Hz')
xlim([-1.5 1])
colorbar
axis xy
caxis([0 0.02])
subplot(223)
imagesc(coh60misc5.time,coh60misc5.freq,squeeze(mean(coh60misc5.cohspctrm(soiidx,:,:),1)));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('coherence diode 2 60 Hz')
xlim([-1.5 1])
colorbar
axis xy
caxis([0 0.02])
subplot(224)
imagesc(coh67misc5.time,coh67misc5.freq,squeeze(mean(coh67misc5.cohspctrm(soiidx,:,:),1)));
xlabel('time (s)')
ylabel('frequency (Hz)')
title('coherence diode 2 67 Hz')
xlim([-1.5 1])
colorbar
axis xy
caxis([0 0.02])
print(fig, fullfile(figpth,'coh_soi'),'-dpng')


%% FFT
% select stimulation and basline interval
cfg = [];
cfg.latency = [0.5 1.5];
STIM60 = ft_selectdata(cfg,datapad60);
cfg.latency = [-1.25 -0.25];
BSL60 = ft_selectdata(cfg,datapad60);

% fft
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.channel = {'meggrad', 'misc004', 'misc005'};
cfg.foilim = [50 80];
cfg.tapsmofrq = 2;
FFT60 = ft_freqanalysis(cfg,STIM60);
FFTBSL60 = ft_freqanalysis(cfg,BSL60);
RC60 = FFT60;
RC60.powspctrm = FFT60.powspctrm./FFTBSL60.powspctrm - 1;
RC60.grad = mGrad;
cfg = [];
cfg.method = 'sum';
RC60 = ft_combineplanar(cfg, RC60);

cfg = [];
ft_multiplotER(cfg,RC60)