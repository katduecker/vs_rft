%% Pilot analysis: VS, alpha and RFT

clear all; close all; clc; beep off;
mpth = 'X:\pilot';
megpth = fullfile(mpth,'meg');
rsppth = fullfile(mpth,'responses');
cdpth = 'X:\experiment';
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')
ft_defaults;
subjcode = 'b57a';
load(fullfile(rsppth,[subjcode,'.mat']))
% load trigger
load(fullfile(cdpth,'trigdef.mat'))

load(fullfile(mpth,'results',[subjcode,'_diodecheck.mat']))

% diode check
subplot(211)
plot(FREQ32y60.freq,FREQ32y60.powspctrm(307:308,:))
title('yellow 60 Hz')
subplot(212)
plot(FREQ32t60.freq,FREQ32t60.powspctrm(307:308,:))
title('teal 60 Hz')


% realtive power change
RC32 = FREQ32rft;
RC32.powspctrm = FREQ32rft.powspctrm./FREQ32bsl.powspctrm- 1;
RC16 = FREQ16rft;
RC16.powspctrm = FREQ16rft.powspctrm./FREQ16bsl.powspctrm- 1;
plot(FREQ32rft.freq,FREQ32rft.powspctrm(307:308,:))
% select grad
cfg = [];
cfg.channel = 'MEGGRAD';
FREQ32grad = ft_selectdata(cfg,RC32);
FREQ16grad = ft_selectdata(cfg,RC16);

% combine planar
cfg = [];
cfg.method ='sum';
FREQ32pl = ft_combineplanar(cfg,FREQ32grad);
FREQ16pl = ft_combineplanar(cfg,FREQ16grad);

% plot
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
%cfg.xlim = [55 100]
%FREQ32pl.powspctrm = log(FREQ32pl.powspctrm);
ft_multiplotER(cfg,FREQ32pl)


%% coherence

load(fullfile(mpth,'results',[subjcode,'_coh_setsize.mat']))

save(fullfile('C:\Users\katha\Desktop\PhD\2_alpha_vs\pilot analysis',[subjcode,'_coh_setsize.mat']))

cfg = [];
cfg.channel = 'MEGGRAD';
ERFgrad = ft_selectdata(cfg,ERF16);

cfg = [];
cfg.layout = 'neuromag306all.lay';
% cfg.baseline = [-1.5 0];
ft_multiplotER(cfg,ERFgrad)

cfg = [];
cfg.layout    = 'neuromag306planar.lay';
cfg.baseline    = [-1.5 0];
cfg.baselinetype = 'relchange';
cfg.refchannel  = 'MISC005';
cfg.colorbar = 'yes';
cfg.zlim = [0 2];
ft_multiplotTFR(cfg,TFR32)

cfg = [];
cfg.avgovertime = 'yes';
cfg.channel = 'MEGGRAD';
TFR32avg = ft_selectdata(cfg,TFR32)

cfg = [];
cfg.layout = 'neuromag306planar.lay';
ft_multiplotER(cfg,TFR32avg)


% plot coherence
cfg = [];
cfg.refchannel = 'MISC004';
cfg.parameter = 'cohspctrm';
ft_multiplotTFR(cfg,COH32rft)

%% PLV

load(fullfile(mpth,'results',[subjcode,'_plv_setsize.mat']))

cfg = [];
cfg.layout = 'neuromag306all.lay';
cfg.parameter = 'plvspctrm';
cfg.refchannel = 'MISC004';
ft_multiplotER(cfg,PLV32rft)