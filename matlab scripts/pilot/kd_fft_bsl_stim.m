function [BSL, STIM, fftbsl, fftstim] = kd_fft_bsl_stim(data,soi,diode1,diode2,bsltoi,stimtoi,avg,foi,gradstruct)

% baseline
cfg = [];
cfg.channel = [soi, diode1,diode2];
cfg.latency = bsltoi;
cfg.avgoverrpt = avg;
BSL = ft_selectdata(cfg,data);

% stim
cfg.latency = stimtoi;
STIM = ft_selectdata(cfg,data);

% FFT
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = foi;
%cfg.pad = 'nextpow2';
fftbsl = ft_freqanalysis(cfg,BSL);
fftstim = ft_freqanalysis(cfg,STIM);
fftbsl.grad = gradstruct;
fftstim.grad = gradstruct;
% combine planar
cfg = [];
cfg.method = 'sum';
fftbsl = ft_combineplanar(cfg,fftbsl);
fftstim = ft_combineplanar(cfg,fftstim);