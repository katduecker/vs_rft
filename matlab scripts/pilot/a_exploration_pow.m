%% Pilot analysis: VS, alpha and RFT

%% Settings
clear all; close all; clc; beep off;
mpth = 'W:\Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
addpath('X:\fieldtrip')
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



%% Separate into 60 and 67 Hz stimulation

% grad structure
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
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 52;
    cfg.detrend = 'yes';
    
    % target tagged at 60 Hz
    cfg.trl = [trl60_misc4;trl60_misc5];  
    block60{fl} = ft_preprocessing(cfg);
    
     % target tagged at 60 Hz
    cfg.trl = [trl67_misc4;trl67_misc5];  
    block67{fl} = ft_preprocessing(cfg);  
end

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
clear data60 data67

trigchan = 'STI101';
misc4 = 'MISC004';
misc5 = 'MISC005';
fig = figure;
plot(ERF60.time,ERF60.avg(strcmp(ERF60.label,misc4),:))
hold on
plot(ERF67.time,ERF67.avg(strcmp(ERF67.label,misc4),:))
xlim([-.2 0.2])
legend('diode 1 60 Hz', 'diode1 67 Hz')
print(fig, fullfile(figpth,'ERF_diode'),'-dpng')
close all

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
cfg.toi = -1.5:0.05:4;
cfg.keeptrials = 'no';
TFR60ind = ft_freqanalysis(cfg,datapad60);
TFR67ind = ft_freqanalysis(cfg,datapad67);
TFR60evo = ft_freqanalysis(cfg,ERF60);
TFR67evo = ft_freqanalysis(cfg,ERF67);
TFR60ind.grad = mGrad;
TFR67ind.grad = mGrad;
TFR60evo.grad = mGrad;
TFR67evo.grad = mGrad;

fig = figure;
cfg = [];
cfg.channel = 'MISC004';

subplot(221)
ft_singleplotTFR(cfg,TFR60ind);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('diode induced 60 Hz')
subplot(222)
ft_singleplotTFR(cfg,TFR67ind);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('diode induced 67 Hz')

subplot(223)
ft_singleplotTFR(cfg,TFR60evo);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('diode evoked 60 Hz')

subplot(224)
ft_singleplotTFR(cfg,TFR67evo);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('diode evoked 67 Hz')
print(fig, fullfile(figpth,'TFR_misc4'),'-dpng')


% combine planar
cfg = [];
cfg.method = 'sum';
TFR60ind = ft_combineplanar(cfg,TFR60ind);
TFR67ind = ft_combineplanar(cfg,TFR67ind);
TFR60evo = ft_combineplanar(cfg,TFR60evo);
TFR67evo = ft_combineplanar(cfg,TFR67evo);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [40 80];
cfg.baseline = [-1.5 -0.5];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
%cfg.xlim = [0.5 4];
fig = figure;
ft_multiplotTFR(cfg,TFR60ind)
print(fig, fullfile(figpth,'TFR60_multi'),'-dpng')
close all
cfg.channel = {'MEG2032+2033', 'MEG2042+2043', 'MEG2112+2113'};
fig = figure;
subplot(221)
ft_singleplotTFR(cfg,TFR60ind);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('60 Hz induced')
subplot(222)
ft_singleplotTFR(cfg,TFR67ind);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('67 Hz induced')
subplot(223)
ft_singleplotTFR(cfg,TFR60evo);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('60 Hz evoked')
caxis([-5 5])
subplot(224)
ft_singleplotTFR(cfg,TFR67evo);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('67 Hz evoked')
caxis([-5 5])
print(fig, fullfile(figpth,'TFR_soi'),'-dpng')


%cfg.channel = 'meggrad';
fig = figure;
cfg = [];
cfg.channel = TFR60ind.label(1:102);
TFR60ind = ft_selectdata(cfg,TFR60ind);
TFR60evo = ft_selectdata(cfg,TFR60evo);
TFR67ind = ft_selectdata(cfg,TFR67ind);
TFR67evo = ft_selectdata(cfg,TFR67evo);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [40 80];
cfg.baseline = [-1.5 -0.5];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
subplot(221)
ft_singleplotTFR(cfg,TFR60ind);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('60 Hz induced')
%caxis([-1 1])
subplot(222)
ft_singleplotTFR(cfg,TFR67ind);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('67 Hz induced')
%caxis([-1 1])
subplot(223)
ft_singleplotTFR(cfg,TFR60evo);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('60 Hz evoked')
%caxis([-5 5])
subplot(224)
ft_singleplotTFR(cfg,TFR67evo);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('67 Hz evoked')
%caxis([-5 5])
print(fig, fullfile(figpth,'TFR_allgrad'),'-dpng')

%% FFT
% select stimulation and basline interval
cfg = [];
cfg.latency = [0.5 1.5];
STIM60 = ft_selectdata(cfg,datapad60);
STIMERF60 = ft_selectdata(cfg,ERF60); 
cfg.latency = [-1.25 -0.25];
BSL60 = ft_selectdata(cfg,datapad60);
BSLERF60 = ft_selectdata(cfg,ERF60); 

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



%%%% OLD
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

% load trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    trl = kd_trlfun(cfg);
    cfg.hpfilter = 'yes';
    cfg.hpfreq   = 52;
    cfg.detrend = 'yes';
    
    cfg.trl = trl;  
    block{fl} = ft_preprocessing(cfg);

end

data = ft_appenddata([],block{:});

% length data
l = cell2mat(cellfun(@length,data.time,'UniformOutput',false));

% discard strange trials
winl = 1;
cfg = [];
cfg.trials = l./fs < 10;
data = ft_selectdata(cfg,data);


%% pad to 8 seconds
datapad = data;
% change time vector
datapad.time(:) = {-2.5:1/fs:5.5-1/fs};
padlength = 8*fs;
% apply zero padding
for t = 1:length(data.trial)
    pad = zeros(length(data.label),padlength-size(data.trial{t},2));
    datapad.trial{t} = horzcat(datapad.trial{t},pad);
end

% plot trigger channel & diode
cfg = [];
ERF = ft_timelockanalysis(cfg,datapad);

trigchan = 'STI101';
misc4 = 'MISC004';
misc5 = 'MISC005';
plot(ERF.time,ERF.avg(strcmp(ERF.label,misc5),:))

xlim([-1.5 2])
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
cfg.toi = -1.5:0.05:4;
cfg.keeptrials = 'no';
TFR = ft_freqanalysis(cfg,ERF);


cfg = [];
cfg.layout = 'neuromag306planar.lay';
cfg.ylim = [40 80];
cfg.baseline = [-1.5 -0.5];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.xlim = [0.5 4];
ft_multiplotTFR(cfg,TFR)