%% Sanity check scripts Visual Search
% explore coherence

clear all; close all; clc; beep off;

% settings
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'data');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts');
figpth = fullfile(mpth,'results','plots','sanity');
pthout = fullfile(mpth,'results','meg','0 sanity');
mkdir(figpth)
addpath(genpath(mtlpth))
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
%addpath('Z:\fieldtrip')
ft_defaults;
subjcode = 'b3ec';
soi = 'MEG2122';

d = dir(megpth);
folds = {d.name};
folds = folds(strncmp(folds,'202',2));
f = cell2mat(cellfun(@(x) strcmp(x(end-3:end),subjcode),folds,'UniformOutput',false));
subpth = fullfile(megpth,folds{f},'meg');
d = dir(subpth);
files = {d.name};
files(1:2) = [];
files = files(cell2mat(cellfun(@(x) strcmp(x(end-6:end),'sss.fif'),files,'UniformOutput',false)));

clear d f folds

% average grad positions

% load in grad structures
grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(subpth,files{fl}))];
end
% average
mGrad = grad(1);
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
clear grad

%% Load in data

% load trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(subpth,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    % separate trials based on target tagging
    cfg.tag = 't';
    cfg.detrend = 'yes';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 54;
    cfg.channel =  {'MEGGRAD','MISC004','MISC005'};
    % load in all trials
    [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg); 
    cfg.trl = trl60_misc4;
    block60misc4{fl} = ft_preprocessing(cfg);
    cfg.trl = trl60_misc5;   
    block60misc5{fl} = ft_preprocessing(cfg);
    cfg.trl = trl67_misc4;
    block67misc4{fl} = ft_preprocessing(cfg);
    cfg.trl = trl67_misc5;
    block67misc5{fl} = ft_preprocessing(cfg);
end

data60misc4 = ft_appenddata([],block60misc4{:});
data60misc4.fsample = block60misc4{1}.fsample;

data60misc5 = ft_appenddata([],block60misc5{:});
data60misc5.fsample = block60misc4{1}.fsample;

data67misc4 = ft_appenddata([],block67misc4{:});
data67misc4.fsample = block60misc4{1}.fsample;

data67misc5 = ft_appenddata([],block67misc5{:});
data67misc5.fsample = block60misc4{1}.fsample;

clear block*

data = ft_appenddata([],data60misc4,data60misc5,data67misc4,data67misc5);
fs = data60misc4.fsample;
% identify reaction time (-3.5 s -> baseline + padding on either side)
rt = (cell2mat(cellfun(@length,data.time,'UniformOutput',false))-4.5*fs)./fs;

% cut out -1 to mean rt + 1
cfg = [];
cfg.latency = [-2 round(mean(rt),2)+1];
data60misc4 = ft_selectdata(cfg,data60misc4);
data60misc5 = ft_selectdata(cfg,data60misc5);
data67misc4 = ft_selectdata(cfg,data67misc4);
data67misc5 = ft_selectdata(cfg,data67misc5);

% select target diodes and concatenate data
cfg = [];
cfg.channel = 'MISC004';
diode60misc4 = ft_selectdata(cfg, data60misc4);
cfg.channel = 'MISC005';
diode60misc5 = ft_selectdata(cfg,data60misc5);
diode60misc5.label = {'diode'};
diode60misc4.label = {'diode'};

% Target tagged at 60 Hz
diode60 = ft_appenddata([],diode60misc4,diode60misc5);
cfg.channel = 'MEGGRAD';
data60misc4 = ft_selectdata(cfg,data60misc4);
data60misc5 = ft_selectdata(cfg,data60misc5);

data60 = ft_appenddata([],data60misc4,data60misc5);
% repeat for 67

% find soi



% % zero pad before averaging
% fs = data60misc4.fsample;
% t = -2.5:1/fs:5.5-1/fs;
% data60misc4 = kd_datapad(data60misc4,fs,8,t);
% data60misc5 = kd_datapad(data60misc5,fs,8,t);
% data67misc4 = kd_datapad(data67misc4,fs,8,t);
% data67misc5 = kd_datapad(data67misc5,fs,8,t);


%% Coherence fieldtrip functions

% 60 Hz tagging
[coh60, psdmeg60, psdmisc60, csd60] = kd_coh_hilb(data60,diode60, foi, frqwdth);
save(fullfile(pthout,[subjcode,'_coh60_mtlb_freqw_',num2str(frqwdth),'_ft_no0pad_sss_T60.mat']),'coh60', 'psdmeg60', 'psdmisc60','csd60')

timevec = data60.time{1};
% PSD
fig = figure;
[~,soiidx] = intersect(data60.label,soi);
subplot(211)
imagesc(timevec,foi,squeeze(psdmisc60(1,:,:)))
title('PSD MISC')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
colorbar
subplot(212)
imagesc(timevec,foi,squeeze(psdmeg60(soiidx,:,:)))
title('PSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
colorbar
print(fig, fullfile(figpth,['psd_meg_diode_60T_ft_frqw_',num2str(frqwdth),'_',soi,'_nohp_no0pad_sss_T60']),'-dpng')

% coherence & CSD% plot coherence and CSD
fig = figure;
subplot(211)
imagesc(timevec,foi,squeeze(csd60(1,:,:)))
title('CSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
%xlim([-1.5 2])
colorbar
subplot(212)
imagesc(timevec,foi,squeeze(coh60(soiidx,:,:)))
title('coherence 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
%xlim([-1.5 2])
colorbar
print(fig, fullfile(figpth,['csd_coh_ft_frqw_',num2str(frqwdth),'_',soi,'_no_hp_no0pad_sss_T60']),'-dpng')



% 60 Hz tagging per diode
frqwdth = 2;                  % width passband     
foi = 55:73;
%[coh60misc4, psdmeg60misc4, psdmisc60misc4, csd60misc4] = kd_coh_hilb(data60misc4, foi, frqwdth,'MISC004', 8);


save(fullfile(pthout,[subjcode,'_coh_mtlb_freqw_',num2str(frqwdth),'_ft_no0pad_sss.mat']),'coh60misc4', 'psdmeg60misc4', 'psdmisc60misc4','csd60misc4')

timevec = data60misc4.time{1};
% PSD
fig = figure;
[~,soiidx] = intersect(data60misc4.label,soi);
subplot(211)
imagesc(timevec,foi,squeeze(psdmisc60misc4(1,:,:)))
title('PSD MISC')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
colorbar
subplot(212)
imagesc(timevec,foi,squeeze(psdmeg60misc4(soiidx,:,:)))
title('PSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
colorbar
print(fig, fullfile(figpth,['psd_meg_diode_ft_frqw_',num2str(frqwdth),'_',soi,'_nohp_no0pad_sss']),'-dpng')

% coherence & CSD% plot coherence and CSD
fig = figure;
[~,soiidx] = intersect(data60misc4.label,soi);
subplot(211)
imagesc(timevec,foi,squeeze(csd60misc4(1,:,:)))
title('CSD MEG SOI')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
%xlim([-1.5 2])
colorbar
subplot(212)
imagesc(timevec,foi,squeeze(coh60misc4(soiidx,:,:)))
title('coherence 60 Hz')
xlabel('time (s)')
ylabel('frequency (Hz)')
axis xy
%xlim([-1.5 2])
colorbar
print(fig, fullfile(figpth,['csd_coh_ft_frqw_',num2str(frqwdth),'_',soi,'_no_hp_no0pad_sss']),'-dpng')


% % save data struct template
% datastruct = data60misc4;
% datastruct = rmfield(datastruct, 'trial');
% datastruct = rmfield(datastruct,'time');
% save(fullfile(mtlpth,'templ_datastruct.mat'),'datastruct','soi')