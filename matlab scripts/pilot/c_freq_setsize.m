%% Pilot analysis: VS, alpha and RFT
% pilot 4: eyes fixated
% SNR set size

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
% load trials
for fl = 1:2%:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    cfg.setsize = 16;
    trl16 = kd_trlfun_trl_setsize(cfg);
    cfg.detrend = 'yes';
    % set size 16
    cfg.trl = trl16;  
    block16{fl} = ft_preprocessing(cfg);
    cfg.setsize = 32;
    trl32 = kd_trlfun_trl_setsize(cfg);
    cfg.trl = trl32;  
    block32{fl} = ft_preprocessing(cfg);
end
block16p = ft_appenddata([],block16{:});
block32p = ft_appenddata([],block32{:});
% find reaction time
rt16 = (cell2mat(cellfun(@length,block16p.time,'UniformOutput',false)) - 2.5*fs - fs)./fs;
rt32 = (cell2mat(cellfun(@length,block32p.time,'UniformOutput',false))- 2.5*fs - fs)./fs;
mrt16 = mean(rt16)
medrt16 = median(rt16)
mrt32 = mean(rt32)
medrt32 = median(rt32)

% change time vector
block16p.time(:) = {-2.5:1/fs:5.5-1/fs};
block32p.time(:) = {-2.5:1/fs:5.5-1/fs};
padlength = 8*fs;
% apply zero padding
% target tagged at 60
for t = 1:length(block16p.trial)
    pad = zeros(length(block16p.label),padlength-size(block16p.trial{t},2));
    block16p.trial{t} = horzcat(block16p.trial{t},pad);
end

% apply zero padding
% target tagged at 60
for t = 1:length(block32p.trial)
    pad = zeros(length(block32p.label),padlength-size(block32p.trial{t},2));
    block32p.trial{t} = horzcat(block32p.trial{t},pad);
end

% average over time
cfg = [];
ERF16 = ft_timelockanalysis(cfg,block16p);
ERF32 = ft_timelockanalysis(cfg,block32p);

% select baseline
cfg = [];
cfg.latency = [-1.25 -0.25];
BSL16 = ft_selectdata(cfg,ERF16);
cfg.avgoverrpt = 'no';
BSL32 = ft_selectdata(cfg,ERF32);

% select flicker
cfg.latency = [0.25 1.25];
FLICK16 = ft_selectdata(cfg,ERF16);
FLICK32 = ft_selectdata(cfg,ERF32);


winl = 1;

%FFT
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.channel = {'meggrad','MISC004','MISC005'};
cfg.foilim = [4 80];
cfg.tapsmofrq = 2;
FFTBSL16 = ft_freqanalysis(cfg,BSL16);
FFTBSL32 = ft_freqanalysis(cfg,BSL32);
FFT16 = ft_freqanalysis(cfg,FLICK16);
FFT32 = ft_freqanalysis(cfg,FLICK32);
FFTBSL16.grad = mGrad;
FFTBSL32.grad = mGrad;
FFT16.grad = mGrad;
FFT32.grad = mGrad;

% diodes% diodes
cfg = [];
cfg.channel = 'MISC004';
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFT16);
title('diode 1 set size 16')
subplot(222)
ft_singleplotER(cfg,FFT16);
title('diode 2 set size 16')
cfg.channel = 'MISC005';
subplot(223)
ft_singleplotER(cfg,FFT32);
title('diode 1 set size 32')
subplot(224)
ft_singleplotER(cfg,FFT32);
title('diode 2 set size 32')
print(fig,fullfile(figpth,'fft_diode_setsize'),'-dpng')

% combine planar
cfg = [];
cfg.method = 'sum';
FFTBSL16cmb = ft_combineplanar(cfg,FFTBSL16);
FFTBSL32cmb = ft_combineplanar(cfg,FFTBSL32);
FFT16cmb = ft_combineplanar(cfg,FFT16);
FFT32cmb = ft_combineplanar(cfg,FFT32);

% relative change
RC16 = FFT16cmb;
RC16.powspctrm = FFT16cmb.powspctrm./FFTBSL16cmb.powspctrm - 1;

RC32 = FFT32cmb;
RC32.powspctrm = FFT32cmb.powspctrm./FFTBSL32cmb.powspctrm - 1;

% plot
% tagging signal
cfg = [];
cfg.channel = {'MEG2032+2033', 'MEG2042+2043'};
cfg.xlim = [50 80];
fig = figure;
subplot(221)
ft_singleplotER(cfg,FFT16cmb)
title('FFT set size 16 evoked')
ylabel('power (T/M)²')
subplot(222)
ft_singleplotER(cfg,RC16)
title('Relative change 16')
ylabel('relative power change')
subplot(223)
ft_singleplotER(cfg,FFT32cmb)
title('FFT set size 16')
ylabel('power (T/M)²')
subplot(224)
ft_singleplotER(cfg,RC32)
title('Relative change 16')
ylabel('relative power change')

print(fig,fullfile(figpth,'fft_soi_rft_setsize'),'-dpng')
