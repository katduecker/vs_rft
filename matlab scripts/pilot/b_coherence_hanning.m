%% Pilot analysis: coherence using complex hanning

% coherence
%% Settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts','pilot');
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

% coherence
%% Settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts','pilot');
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
% tagging T and D
tagblock = {'t','t'};
% load trials
for fl = 1:2%:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    
    % separate into trials based on the tagging of the target
    cfg.tag = tagblock{fl};
    [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg);
    cfg.bpfilter = 'yes';
    cfg.bpfreq   = [52 75];
    cfg.detrend = 'yes';
    cfg.channel = {'MEGGRAD','MISC004','MISC005'};
    
    % misc 4 tagged at 60 Hz, misc5 at 67 Hz
    cfg.trl = [trl60_misc4; trl67_misc5];
    misc4_60_misc5_67{fl} = ft_preprocessing(cfg);
    % tagging at 60 Hz
    cfg.trl = [trl60_misc5; trl67_misc4];
    misc5_60_misc4_67{fl} = ft_preprocessing(cfg);
end

data60misc4 = ft_appenddata([],misc4_60_misc5_67{:});
data60misc5 = ft_appenddata([],misc5_60_misc4_67{:});
data60misc4.fsample = 1000;
data60misc5.fsample = 1000;

% Coherence
foi = 55:62;
N = 3;
[coh_spct, psd_meg, psd_misc,csd_meg_misc] = kd_coh_hann(data60misc4, foi, 'MISC004', data60misc4.fsample, N, 8);

soi = 'MEG2043';
[~,soiidx] = intersect(data60misc4.label,soi);
fig = figure;
imagesc(linspace(-2.5,5.5,8000),56:62,squeeze(psd_meg(soiidx,:,:)))
title('PSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
xlim([-1.5 4])
colorbar
print(fig, fullfile(figpth,['coh_hann_N_',num2str(N)]),'-dpng')
