%% Pilot analysis: VS, alpha and RFT

% plot coherence
clear all; close all; clc; beep off;
mpth = 'X:\';
cd(mpth)
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')
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
soi = {'MEG2032','MEG2033','MEG2112','MEG2113'};


load(fullfile(pthout,subjcode,[subjcode,'coh_hilbert_ftrip_fwidth_3.mat']))

% multiplot
cfg = [];
cfg.parameter = 'cohspctrm';
%cfg.channel = 'MEGGRAD';
cfg.refchannel = 'MISC004';
cfg.layout = 'neuromag306planar.lay';
cfg.colorbar = 'yes';
cfg.ylim = [50 80];
ft_multiplotTFR(cfg,coh_meg_misc4)
ft_hastoolbox('brewermap',1);
colormap(flipud(brewermap(64,'RdBu')))

for s = 1:length(soi)
    subplot(2,2,s)
    idx = find(strcmp(coh_meg_misc4.label,soi{s}));
    imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(coh_meg_misc4.cohspctrm(idx,:,:)))
    axis xy
    xlim([-1.5 1])
    ylim([55 75])
    cb = colorbar;
    cb.Label.String = 'coherence';
    title(soi{s})
    xlabel('time (s)')
    ylabel('freq (Hz)')
end

print(fullfile(figpth,'coh_soi'),'-dpng')
close all
% sanity: plot psd of misc

% here it doesn't matter because all rows in psd_misc4 are the same..
misc4 = find(strcmp(coh_meg_misc4.label,'MISC004'));
misc5 = find(strcmp(coh_meg_misc4.label,'MISC005'));
subplot(211)
imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(psd_misc4(misc4,:,:)))
axis xy
xlim([-1.5 1])
xlabel('time (s)')
ylabel('freq (Hz)')
title('MISC004')
cb = colorbar;
cb.Label.String = 'power';
subplot(212)
imagesc(coh_meg_misc4.time,coh_meg_misc4.freq,squeeze(psd_meg(misc5,:,:)))
axis xy
xlim([-1.5 1])
xlabel('time (s)')
ylabel('freq (Hz)')
cb = colorbar;
cb.Label.String = 'power';
title('MISC005')
print(fullfile(figpth,'psd_msic'),'-dpng')

tmpl = coh_meg_misc4;           % use as template
clear coh* csd* psd*
%% 
subjcode = 'b59c';
load(fullfile(pthout,subjcode,[subjcode,'coh_hilbert_ftrip_fwidth_3.mat']))
figpth = fullfile(mpth,'pilot','results','plots',subjcode);
mkdir(figpth)
tmpl.time = tmpl.time{1};
tmpl.time = tmpl.time(1:length(tmpl.time)-1);

subplot(211)
imagesc(tmpl.time,tmpl.freq,squeeze(psd_misc4(misc4,:,:)))
axis xy
xlim([-1.5 1])
xlabel('time (s)')
ylabel('freq (Hz)')
title('surrog 60 Hz')
cb = colorbar;
cb.Label.String = 'power';
subplot(212)
imagesc(tmpl.time,tmpl.freq,squeeze(psd_meg(misc5,:,:)))
axis xy
xlim([-1.5 1])
xlabel('time (s)')
ylabel('freq (Hz)')
cb = colorbar;
cb.Label.String = 'power';
title('surrog 67 Hz')
print(fullfile(figpth,'psd_msic'),'-dpng')


% coherence in template
tmpl.cohspctrm = coh_meg_misc4_spct;

% multiplot
cfg = [];
cfg.parameter = 'cohspctrm';
%cfg.channel = 'MEGGRAD';
cfg.refchannel = 'MISC004';
cfg.layout = 'neuromag306planar.lay';
cfg.colorbar = 'yes';
cfg.ylim = [50 80];
ft_multiplotTFR(cfg,tmpl)
ft_hastoolbox('brewermap',1);
colormap(flipud(brewermap(64,'RdBu')))


for s = 1:length(soi)
    subplot(2,2,s)
    idx = find(strcmp(tmpl.label,soi{s}));
    imagesc(tmpl.time,tmpl.freq,squeeze(tmpl.cohspctrm(idx,:,:)))
    axis xy
    xlim([-1.5 1])
    ylim([55 75])
    cb = colorbar;
    cb.Label.String = 'coherence';
    title(soi{s})
    xlabel('time (s)')
    ylabel('freq (Hz)')
end
print(fullfile(figpth,'coh_soi'),'-dpng')
