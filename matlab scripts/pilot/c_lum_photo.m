% Equipment test: luminance of colours

% settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
ft_defaults;

% path to fif files
pathin = fullfile(mpth,'luminance_diode','201207');

% colours
clrs = {'pink','teal','yellow'};

misc4 = 'MISC004';

% load in fif files
for c = 1:length(clrs)
    cfg = [];
    cfg.dataset  = fullfile(pathin,[clrs{c},'.fif']);
    % get trial specifics
    cfg.path = cdpth;               % experiment path with trigger info
    [~, ~, trl4,~, ~, trl5] = kd_trlfun_trl_phd(cfg);
    %trl = kd_trlfun(cfg);
    cfg.detrend = 'yes';
    cfg.trl = trl4(2,:);
    cfg.channel = misc4;
    phtd60{c} = ft_preprocessing(cfg);          % misc4 picks up 60 Hz signal
        
    cfg.trl = trl5(2,:);
    phtd67{c} = ft_preprocessing(cfg);          % misc4 picks up 67 Hz signal
end

% cut out stimulation time
cfg = [];
cfg.latency = [-.5 4];
for c = 1:length(clrs)
    phtd60{c} = ft_selectdata(cfg,phtd60{c});
    phtd67{c} = ft_selectdata(cfg,phtd67{c});
end

% check timing
for c = 1:length(clrs)
    for t = 1:length(phtd60{c}.trial)
        plot(phtd60{c}.time{t},phtd60{c}.trial{t})
        hold on
        xlim([-0.5 0.5])
    end
end

% frequency analysis
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.taper = 'hanning';
cfg.foilim = [57 70];
cfg.tapsomfrq = 0;
for c = 1:length(clrs)
    fft60{c} = ft_freqanalysis(cfg,phtd60{c});
    fft67{c} = ft_freqanalysis(cfg,phtd67{c});  
end


plot(fft60{1}.freq,fft60{1}.powspctrm,'Color','m')
hold on
plot(fft60{1}.freq,fft60{2}.powspctrm,'Color','c')
hold on
plot(fft60{1}.freq,fft60{3}.powspctrm,'Color','y')

hold on
plot(fft67{1}.freq,fft67{1}.powspctrm,'Color','m')
hold on
plot(fft67{1}.freq,fft67{2}.powspctrm,'Color','c')
hold on
plot(fft67{1}.freq,fft67{3}.powspctrm,'Color','y')
xlabel('frequency (Hz)')
ylabel('luminance (photodiode)')
%title('luminance of selected fully saturated colours')

[~,p60] = min(abs(fft60{1}.freq-60));
[~,p67] = min(abs(fft60{1}.freq-67));

ampl = [];
for c = 1:length(clrs)
    ampl = [ampl; fft60{c}.powspctrm(p60),fft67{c}.powspctrm(p67)];
end