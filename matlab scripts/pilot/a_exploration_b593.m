%% Pilot analysis: VS, alpha and RFT
% pilot 4: eyes fixated

% frequency analyses
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
subjcode = 'b593';
winl = 1;
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

%% Use this trial function as the others are not working..
% tagging T and D
tagblock = {'t','t'};
% load trials
for fl = 1:2%:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    
    % target or distractor tagging?
    cfg.tag = tagblock{fl};
    [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg);
%     cfg.hpfilter = 'yes';
%     cfg.hpfreq   = 52;
    cfg.detrend = 'yes';
    
    % tagging at 60 Hz
    cfg.trl = [trl60_misc4;trl60_misc5;trl67_misc4;trl67_misc5];  
    block{fl} = ft_preprocessing(cfg);
end

data = ft_appenddata([],block{:});

% identify reaction time (-4.5 s -> baseline + padding on either side)
rt = (cell2mat(cellfun(@length,data.time,'UniformOutput',false))-4.5*fs)./fs;

% select trials with rt > 0.75 (to have at least 250 ms after ERF)
cfg = [];
cfg.trials = logical(rt > 0.75);
data = ft_selectdata(cfg,data);

rt = rt(rt > 0.75);
mrt = mean(rt);

% cut out data up until mean rt
cfg = [];
cfg.latency = [-1.5 mrt+winl/2];
data = ft_selectdata(cfg,data);

% datapad = data;
% % change time vector
% datapad.time(:) = {-2.5:1/fs:5.5-1/fs};
% padlength = 8*fs;
% % apply zero padding
% % target tagged at 60
% for t = 1:length(datapad.trial)
%     pad = zeros(length(data.label),padlength-size(data.trial{t},2));
%     datapad.trial{t} = horzcat(datapad.trial{t},pad);
% end

% plot trigger channel & diode
cfg = [];
ERF = ft_timelockanalysis(cfg,data);
% trigchan = 'STI101';
% misc4 = 'MISC004';
% misc5 = 'MISC005';
% plot(ERF.time,ERF.avg(strcmp(ERF.label,trigchan),:))
% xlim([-2 1])

% frequency analysis
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
% cfg.pad = 8;
% cfg.padtype = 'mean';
cfg.foi = 4:1:80;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
cfg.toi = -1.25:0.05:1.25;
cfg.keeptrials = 'no';
TFR_ev = ft_freqanalysis(cfg,ERF);
%TFR_ind = ft_freqanalysis(cfg,data);
TFR_ev.grad = mGrad;
%TFR_ind.grad = mGrad;
% 
% % diodes
% fig = figure;
% cfg = [];
% cfg.channel = 'MISC004';
% subplot(221)
% ft_singleplotTFR(cfg,TFR_ev);
% title('diode 1 evoked')
% subplot(222)
% ft_singleplotTFR(cfg,TFR_ind);
% title('diode 1 induced')
% cfg.channel = 'MISC005';
% subplot(223)
% ft_singleplotTFR(cfg,TFR_ev);
% title('diode 2 evoked')
% subplot(224)
% ft_singleplotTFR(cfg,TFR_ind);
% title('diode 2 induced')
% print(fig,fullfile(figpth,'tfr_diode'),'-dpng')
% close all

cfg = [];
cfg.method = 'sum';
%TFR_indcmb = ft_combineplanar(cfg,TFR_ind);
TFR_evcmb = ft_combineplanar(cfg,TFR_ev);

% cfg = [];
% cfg.layout = 'neuromag306cmb.lay';
% cfg.ylim = [55 70];
% cfg.baseline = [-1.5+winl/2 0-winl/2];
% cfg.baselinetype = 'relchange';
% cfg.colorbar = 'yes';
% %cfg.zlim = [0 8];
% %cfg.xlim = [0.5 4];
% ft_multiplotTFR(cfg,TFR_evcmb)

soi = 'MEG2032+2033';

% average TFR -> spectrum
cfg = [];
cfg.baseline =  [-1.5+winl/2 0-winl/2];
cfg.baselinetype = 'relchange';
TFR_evcmb_rc = ft_freqbaseline(cfg,TFR_evcmb);

cfg = [];
cfg.latency = [0.5 mrt];
cfg.avgovertime = 'yes';
cfg.channel = soi;
AVG_TFRev = ft_selectdata(cfg,TFR_evcmb_rc);

% plot
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [55 70];
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.channel = soi;
cfg.zlim = [-8 8];
cfg.xlim = [-1 1];
fig = figure;
subplot(211)
ft_singleplotTFR(cfg,TFR_evcmb)
title(['TFR induced ', soi])
xlabel('time (s)')
ylabel('frequency (Hz)')
cb = colorbar;
cb.Label.String = 'relative power change';
subplot(212)
%cfg.zlim = [0 8];
plot(AVG_TFRev.freq,AVG_TFRev.powspctrm);
xlim([55 75])
title(['averaged TFR evoked, ',soi])
xlabel('time (s)')
ylabel('frequency (Hz)')

print(fig,fullfile(figpth,['tfr_',soi,'_winl',num2str(winl*1000)]),'-dpng')
close all
%% FFT

% select baseline
cfg = [];
cfg.latency = [-1.25 -0.25];
cfg.avgoverrpt = 'yes';
BSLev = ft_selectdata(cfg,datapad);
cfg.avgoverrpt = 'no';
BSLind = ft_selectdata(cfg,datapad);

% select flicker
cfg.latency = [0.25 1.25];
FLICKind = ft_selectdata(cfg,datapad);
cfg.avgoverrpt = 'yes';
FLICKev = ft_selectdata(cfg,datapad);

% 
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.channel = {'meggrad','MISC004','MISC005'};
cfg.foilim = [4 80];
cfg.tapsmofrq = 2;
FFTbslev = ft_freqanalysis(cfg,BSLev);
FFTflickev = ft_freqanalysis(cfg,FLICKev);
FFTbslind = ft_freqanalysis(cfg,BSLind);
FFTflickind = ft_freqanalysis(cfg,FLICKind);
FFTbslev.grad = mGrad;
FFTflickev.grad = mGrad;
FFTbslind.grad = mGrad;
FFTflickind.grad = mGrad;

% diodes
cfg = [];
cfg.channel = 'MISC004';
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFTflickev);
title('diode 1 evoked')
subplot(222)
ft_singleplotER(cfg,FFTflickind);
title('diode 1 induced')
cfg.channel = 'MISC005';
subplot(223)
ft_singleplotER(cfg,FFTflickev);
title('diode 2 evoked')
subplot(224)
ft_singleplotER(cfg,FFTflickind);
title('diode 2 induced')
print(fig,fullfile(figpth,'fft_diode'),'-dpng')
close all

% combine planar
cfg = [];
cfg.method = 'sum';
FFTbslevcmb = ft_combineplanar(cfg,FFTbslev);
FFTbslindcmb = ft_combineplanar(cfg,FFTbslind);
FFTflickevcmb = ft_combineplanar(cfg,FFTflickev);
FFTflickindcmb = ft_combineplanar(cfg,FFTflickind);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
ft_multiplotER(cfg,FFTflickevcmb);
ft_multiplotER(cfg,FFTflickindcmb);

% relative power change
RCflickev = FFTflickevcmb;
RCflickev.powspctrm = FFTflickevcmb.powspctrm./FFTbslevcmb.powspctrm - 1;
RCflickind = FFTflickindcmb;
RCflickind.powspctrm = FFTflickindcmb.powspctrm./FFTbslindcmb.powspctrm - 1;

% tagging signal
cfg = [];
cfg.channel = {'MEG2032+2033', 'MEG2042+2043'};
cfg.xlim = [50 80];
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFTflickevcmb)
title('power spectrum SOI evoked')
ylabel('power (T/M)²')
subplot(222)
ft_singleplotER(cfg,RCflickev)
title('relative power change SOI evoked')
ylabel('relative power change')
subplot(223)
ft_singleplotER(cfg,FFTflickindcmb)
title('power spectrum SOI induced')
ylabel('power (T/M)²')
subplot(224)
ft_singleplotER(cfg,RCflickind)
title('relative power change SOI induced')
print(fig,fullfile(figpth,'fft_soi_rft'),'-dpng')

% alpha
cfg.xlim = [4 30];
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFTflickevcmb)
title('lower frequencies phase-locked')
ylabel('power (T/M)²')
subplot(222)
ft_singleplotER(cfg,RCflickev)
title('relative change phase-locked')
ylabel('relative power change')
subplot(223)
ft_singleplotER(cfg,FFTflickindcmb)
title('low freq not phase-locked')
ylabel('power (T/M)²')
subplot(224)
ft_singleplotER(cfg,RCflickind)
title('relative change not phase-locked')
print(fig,fullfile(figpth,'fft_soi_alpha'),'-dpng')


%% Tagging target only
clear data* FFT* TFR* BSL* FLICK* ERF* RC* trl*

tagblock = {'t','t','t','t'};
% load trials
for fl = 3:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    
    % target or distractor tagging?
    cfg.tag = tagblock{fl};
    [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg);
%     cfg.hpfilter = 'yes';
%     cfg.hpfreq   = 52;
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [52 80];
    cfg.detrend = 'yes';
    
    % tagging at 60 Hz
    cfg.trl = [trl60_misc4;trl60_misc5];  
    block60{fl} = ft_preprocessing(cfg);
    
     % tagging at 67 Hz
    cfg.trl = [trl67_misc4;trl67_misc5];  
    block67{fl} = ft_preprocessing(cfg);  
end
block60(1:2) = [];
block67(1:2) = [];
data60 = ft_appenddata([],block60{:});
data67 = ft_appenddata([],block67{:});

% discard strange trials
l = cell2mat(cellfun(@length,data60.time,'UniformOutput',false));
cfg = [];
cfg.trials = l./fs < 10;
data60 = ft_selectdata(cfg,data60);
% pad data to 8 seconds
datapad60 = data60;
datapad67 = data67;
% change time vector
datapad60.time(:) = {-2.5:1/fs:5.5-1/fs};
datapad67.time(:) = {-2.5:1/fs:5.5-1/fs};
padlength = 8*fs;
% apply zero padding
% target tagged at 60
for t = 1:length(datapad60.trial)
    pad = zeros(length(data60.label),padlength-size(data60.trial{t},2));
    datapad60.trial{t} = horzcat(datapad60.trial{t},pad);
end
% target tagged at 67
for t = 1:length(datapad67.trial)
    pad = zeros(length(data67.label),padlength-size(data67.trial{t},2));
    datapad67.trial{t} = horzcat(datapad67.trial{t},pad);
end

% plot trigger channel & diode
cfg = [];
ERF60 = ft_timelockanalysis(cfg,datapad60);
ERF67 = ft_timelockanalysis(cfg,datapad67);
trigchan = 'STI101';
misc4 = 'MISC004';
misc5 = 'MISC005';
plot(ERF60.time,ERF67.avg(strcmp(ERF67.label,misc5),:))
xlim([-2 1])
% ylim([0 10])



%% TFR
winl = 1;
% frequency analysis
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
% cfg.pad = 8;
% cfg.padtype = 'mean';
cfg.foi = 4:1:80;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
cfg.toi = -1.25:0.05:1.25;
cfg.keeptrials = 'no';
TFR60_ev = ft_freqanalysis(cfg,ERF60);
TFR60_ind = ft_freqanalysis(cfg,datapad60);
TFR60_ev.grad = mGrad;
TFR60_ind.grad = mGrad;
TFR67_ev = ft_freqanalysis(cfg,ERF67);
TFR67_ind = ft_freqanalysis(cfg,datapad67);
TFR67_ev.grad = mGrad;
TFR67_ind.grad = mGrad;

% diodes
fig = figure;
cfg = [];
cfg.channel = 'MISC004';
subplot(221)
ft_singleplotTFR(cfg,TFR60_ev);
title('diode 1 evoked 60')
subplot(222)
ft_singleplotTFR(cfg,TFR60_ind);
title('diode 1 induced 60')
cfg.channel = 'MISC005';
subplot(223)
ft_singleplotTFR(cfg,TFR60_ev);
title('diode 2 evoked 60')
subplot(224)
ft_singleplotTFR(cfg,TFR60_ind);
title('diode 2 induced 60')
print(fig,fullfile(figpth,'tfr_60Hz_ttag'),'-dpng')
close all

% diodes
fig = figure;
cfg = [];
cfg.channel = 'MISC004';
subplot(221)
ft_singleplotTFR(cfg,TFR67_ev);
title('diode 1 evoked 67')
subplot(222)
ft_singleplotTFR(cfg,TFR67_ind);
title('diode 1 induced 67')
cfg.channel = 'MISC005';
subplot(223)
ft_singleplotTFR(cfg,TFR67_ev);
title('diode 2 evoked 67')
subplot(224)
ft_singleplotTFR(cfg,TFR67_ind);
title('diode 2 induced 67')
print(fig,fullfile(figpth,'tfr_67Hz_ttag'),'-dpng')
close all

cfg = [];
cfg.method = 'sum';
TFR60_indcmb = ft_combineplanar(cfg,TFR60_ind);
TFR60_evcmb = ft_combineplanar(cfg,TFR60_ev);
TFR67_indcmb = ft_combineplanar(cfg,TFR67_ind);
TFR67_evcmb = ft_combineplanar(cfg,TFR67_ev);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [55 70];
cfg.baseline = [-1.25 -0.25];
cfg.baselinetype = 'db';
cfg.colorbar = 'yes';
cfg.zlim = 'zeromax';
cfg.channel = {'MEG2032+2033', 'MEG2042+2043'};
fig = figure;
subplot(221)
ft_singleplotTFR(cfg,TFR60_indcmb)
title('TFR induced over SOI 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
cb = colorbar;
cb.Label.String = 'db';
subplot(222)
ft_singleplotTFR(cfg,TFR67_indcmb)
title('TFR induced over SOI 67 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
cb = colorbar;
cb.Label.String = 'db';
subplot(223)
%cfg.zlim = [0 8];
ft_singleplotTFR(cfg,TFR60_evcmb)
title('TFR evoked over SOI 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
cb = colorbar;
cb.Label.String = 'db';
subplot(224)
%cfg.zlim = [0 8];
ft_singleplotTFR(cfg,TFR67_evcmb)
title('TFR evoked over SOI 67 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
cb = colorbar;
cb.Label.String = 'db';
print(fig,fullfile(figpth,'tfr_ttag_soi'),'-dpng')
close all

%% FFT

% select baseline
cfg = [];
cfg.latency = [-1.25 -0.25];
cfg.avgoverrpt = 'yes';
BSLev = ft_selectdata(cfg,datapad60);
cfg.avgoverrpt = 'no';
BSLind = ft_selectdata(cfg,datapad60);

% select flicker
cfg.latency = [0.25 1.25];
FLICKind = ft_selectdata(cfg,datapad60);
cfg.avgoverrpt = 'yes';
FLICKev = ft_selectdata(cfg,datapad60);

% 
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.channel = {'meggrad','MISC004','MISC005'};
cfg.foilim = [4 80];
cfg.tapsmofrq = 2;
FFTbslev = ft_freqanalysis(cfg,BSLev);
FFTflickev = ft_freqanalysis(cfg,FLICKev);
FFTbslind = ft_freqanalysis(cfg,BSLind);
FFTflickind = ft_freqanalysis(cfg,FLICKind);
FFTbslev.grad = mGrad;
FFTflickev.grad = mGrad;
FFTbslind.grad = mGrad;
FFTflickind.grad = mGrad;

% diodes
cfg = [];
cfg.channel = 'MISC004';
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFTflickev);
title('diode 1 evoked 60')
subplot(222)
ft_singleplotER(cfg,FFTflickind);
title('diode 1 induced 60')
cfg.channel = 'MISC005';
subplot(223)
ft_singleplotER(cfg,FFTflickev);
title('diode 2 evoked 60')
subplot(224)
ft_singleplotER(cfg,FFTflickind);
title('diode 2 induced 60')
print(fig,fullfile(figpth,'fft_diode_60ttag'),'-dpng')
close all

% combine planar
cfg = [];
cfg.method = 'sum';
FFTbslevcmb = ft_combineplanar(cfg,FFTbslev);
FFTbslindcmb = ft_combineplanar(cfg,FFTbslind);
FFTflickevcmb = ft_combineplanar(cfg,FFTflickev);
FFTflickindcmb = ft_combineplanar(cfg,FFTflickind);

% relative power change
RCflickev = FFTflickevcmb;
RCflickev.powspctrm = FFTflickevcmb.powspctrm./FFTbslevcmb.powspctrm - 1;
RCflickind = FFTflickindcmb;
RCflickind.powspctrm = FFTflickindcmb.powspctrm./FFTbslindcmb.powspctrm - 1;

% tagging signal
cfg = [];
cfg.channel = {'MEG2032+2033', 'MEG2042+2043'};
cfg.xlim = [50 80];
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFTflickevcmb)
title('power spectrum SOI evoked 60 Hz')
ylabel('power (T/M)²')
subplot(222)
ft_singleplotER(cfg,RCflickev)
title('relative power change SOI evoked')
ylabel('relative power change')
subplot(223)
ft_singleplotER(cfg,FFTflickindcmb)
title('power spectrum SOI induced')
ylabel('power (T/M)²')
subplot(224)
ft_singleplotER(cfg,RCflickind)
title('relative power change SOI induced')
print(fig,fullfile(figpth,'fft_soi_ttag60'),'-dpng')