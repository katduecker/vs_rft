%% Sanity check scripts Visual Search

clear all; close all; clc; beep off;

% settings
mpth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
cd(mpth)
megpth = fullfile(mpth,'data');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab scripts');
figpth = fullfile(mpth,'results','plots','sanity');
mkdir(figpth)
addpath(genpath(mtlpth))
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
%addpath('Z:\fieldtrip')
ft_defaults;
subjcode = 'b3ec';
sens = 'meggrad';
magsoi = {'MEG2111'};
gradsoi = {'MEG2122+2123'};
sss = 0;
winl = 0.5;


d = dir(megpth);
folds = {d.name};
folds = folds(strncmp(folds,'202',2));
f = cell2mat(cellfun(@(x) strcmp(x(end-3:end),subjcode),folds,'UniformOutput',false));
subpth = fullfile(megpth,folds{f},'meg');
d = dir(subpth);
files = {d.name};
files(1:2) = [];
if sss
    files = files(cell2mat(cellfun(@(x) strcmp(x(end-6:end),'sss.fif'),files,'UniformOutput',false)));
else
    files = files(strncmp(files,'part',4));
    files = files(cell2mat(cellfun(@(x) ~strcmp(x(end-6:end),'sss.fif'),files,'UniformOutput',false)));
end
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
    % target or distractor tagging?
    cfg.detrend = 'yes';
    %cfg.hpfilter = 'yes';
    %cfg.hpfreq = 52;
    % load in all trials
    cfg.trl = kd_trlfun_alltrl(cfg); 
    block{fl} = ft_preprocessing(cfg);
end

data = ft_appenddata([],block{:});
fs = block{fl}.fsample;                 % sampling rate
clear block

% find SOI
bsltoi = [-1 -0.5];
stimtoi = [0.25 .75];
tagfreq = 60;

soi60 = kd_fsoi(data,bsltoi,stimtoi,sens,tagfreq,winl,mGrad);

tagfreq = 67;
soi67 = kd_fsoi(data,bsltoi,stimtoi,sens,tagfreq,winl,mGrad);

% ignore for now that they don't overlap perfectly
soicmb = intersect(soi60,soi67);
% pull cmbs apart
soi = {};

for s = 1:length(soicmb)
    soi = [soi,soicmb{s}(1:strfind(soicmb{1},'+')-1),['MEG',soicmb{s}(strfind(soicmb{1},'+')+1:end)]];
end


% pad data to 8 sec
% zero pad before averaging
% t = -2.5:1/fs:5.5-1/fs;
% datapad = kd_datapad(data,fs,8,t);


% plot trigger channel & diode
cfg = [];
cfg.latency = [-1.25 1.25];
cfg.avgoverrpt = 'yes';
ERF = ft_selectdata(cfg,data);

% trigchan = 'STI101';
% misc4 = 'MISC004';
% misc5 = 'MISC005';
% subplot(311)
% for t = 1:length(datapad.trial)
%     plot(datapad.time{1},datapad.trial{t}(strcmp(datapad.label,trigchan),:))
% end
% xlim([-1.75 1])
% ylim([0 50])
% title('trigger')
% subplot(312)
% for t = 1:length(datapad.trial)
%     plot(datapad.time{1},datapad.trial{t}(strcmp(datapad.label,misc4),:))
% end
% xlim([-1.75 1])
% title('diode 1')
% subplot(313)
% for t = 1:length(datapad.trial)
%     plot(datapad.time{1},datapad.trial{t}(strcmp(datapad.label,misc5),:))
% end
% xlim([-1.75 1])

%% Frequency analysis different window length

% frequency analysis
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = soi;
% cfg.pad = 8;
% cfg.padtype = 'mean';
cfg.foi = 56:1/winl:80;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
cfg.toi = -1.25:0.05:1.25;
cfg.keeptrials = 'no';
TFR_ev = ft_freqanalysis(cfg,ERF);  
TFR_ev.grad = mGrad;

if strcmp(sens,'meggrad')
    % find soi
    cfg = [];
    cfg.method = 'sum';
    TFR_ev = ft_combineplanar(cfg,TFR_ev);
end

cfg = [];
if strcmp(sens,'megmag')
    cfg.layout = 'neuromag306mag.lay';
elseif strcmp(sens,'meggrad')
    cfg.layout = 'neuromag306cmb.lay';
end
cfg.ylim = [58 80];
cfg.baseline = [-1.5 -0.5];
cfg.baselinetype ='relchange';
cfg.colorbar = 'yes';
cfg.zlim = [-6 6];
%cfg.xlim = [0.5 4];
ft_multiplotTFR(cfg,TFR_ev)


% average TFR
cfg = [];
cfg.baseline = [-winl/2-1 -winl/2];
cfg.baselinetype = 'relchange';
TFR_evbsl = ft_freqbaseline(cfg,TFR_ev);

cfg = [];
cfg.latency = [0 1];
cfg.avgovertime = 'yes';
AVGev = ft_selectdata(cfg,TFR_evbsl);

fig = figure;
subplot(211)
cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [58 80];
cfg.baseline = [-1.5 -winl/2];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.channel = 8;
if strcmp(sens,'megmag')
    %cfg.channel = magsoi;
    cfg.zlim = [-6 6];
elseif strcmp(sens,'meggrad')
    %cfg.channel = gradsoi;
    %cfg.zlim = [-1.5 1.5];
end
ft_singleplotTFR(cfg,TFR_ev);
xlabel('time (s)')
ylabel('frequency (Hz)')
%title('[all trials soi ',cfg.channel{:},'; evoked TFR')
subplot(212)
%soiidx = ismember(AVGev.label,cfg.channel);
plot(AVGev.freq,AVGev.powspctrm(8,:))
xlim([55 80])
xlabel('frequency (Hz)')
ylabel('power change to baseline')
%title(['all trials soi ',cfg.channel{:},'; averaged evoked TFR 0.5 - 1.5 s'])
if sss
    print(fig,fullfile(figpth,[subjcode,'_averaged_tfr_soi_winl_',num2str(winl*1000),'_sss_',sens]),'-dpng')
else
    print(fig,fullfile(figpth,[subjcode,'_averaged_tfr_soi_winl_',num2str(winl*1000),sens]),'-dpng')
end




% 
% %% Frequency analysis evoked vs induced
% winl = 0.25;
% % frequency analysis
% cfg = [];
% cfg.method = 'mtmconvol';
% cfg.output = 'pow';
% cfg.channel = {'MEGGRAD','MISC004','MISC005'};
% % cfg.pad = 8;
% % cfg.padtype = 'mean';
% cfg.foi = 4:1:80;
% cfg.taper = 'hanning';
% cfg.t_ftimwin = ones(size(cfg.foi)).*winl;
% cfg.toi = -1.25:0.05:1.25;
% cfg.keeptrials = 'no';
% TFR_ev = ft_freqanalysis(cfg,ERF);                  % evoked/phase-locked
% TFR_ind = ft_freqanalysis(cfg,datapad);             % induced/non-phase locked
% TFR_ev.grad = mGrad;
% TFR_ind.grad = mGrad;
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
% print(fig,fullfile(figpth,[subjcode,'_tfr_diode']),'-dpng')
% close all
% 
% cfg = [];
% cfg.method = 'sum';
% TFR_indcmb = ft_combineplanar(cfg,TFR_ind);
% TFR_ev = ft_combineplanar(cfg,TFR_ev);
% 
% cfg = [];
% cfg.layout = 'neuromag306cmb.lay';
% cfg.ylim = [58 80];
% cfg.baseline = [-1.25 -0.1];
% cfg.baselinetype = 'relchange';
% cfg.colorbar = 'yes';
% %cfg.zlim = [0 8];
% %cfg.xlim = [0.5 4];
% ft_multiplotTFR(cfg,TFR_ev)
% 
% % singleplot
% soi = {'MEG2032+2033', 'MEG2042+2043'};
% cfg.channel = soi;
% fig = figure;
% subplot(211)
% ft_singleplotTFR(cfg,TFR_indcmb)
% title('TFR induced over SOI')
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% cb = colorbar;
% cb.Label.String = 'power (T/M)²';
% subplot(212)
% %cfg.zlim = [0 8];
% ft_singleplotTFR(cfg,TFR_ev)
% title('TFR evoked over SOI')
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% cb = colorbar;
% cb.Label.String = 'power (T/M)²';
% print(fig,fullfile(figpth,[subjcode,'_tfr_soi_nobsl_sss']),'-dpng')
% 
% TFR_ind_log = TFR_indcmb;
% TFR_ind_log.powspctrm = log(TFR_ind_log.powspctrm);
% TFR_ev_log = TFR_ev;
% TFR_ev_log.powspctrm = log(TFR_ev_log.powspctrm);
% 
% soi = {'MEG2032+2033', 'MEG2042+2043'};
% cfg.channel = soi;
% cfg.ylim = [52 80];
% fig = figure;
% subplot(211)
% ft_singleplotTFR(cfg,TFR_ind_log)
% title('TFR induced over SOI')
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% cb = colorbar;
% cb.Label.String = 'log(power (T/M)²)';
% subplot(212)
% %cfg.zlim = [0 8];
% ft_singleplotTFR(cfg,TFR_ev_log)
% title('TFR evoked over SOI')
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% cb = colorbar;
% cb.Label.String = 'log(power (T/M)²)';
% if sss
%     print(fig,fullfile(figpth,[subjcode,'_tfr_soi_log_sss']),'-dpng')
% else
%     print(fig,fullfile(figpth,[subjcode,'_tfr_soi_log']),'-dpng')
% end
% 
% 
% %% Power change
% 
% soi = {'MEG2032','MEG2033', 'MEG2042','MEG2043'};
% soicmb = {'MEG2032+2033', 'MEG2042+2043'};
% 
% % evoked
% [~, ~, bslev, rftev] = kd_fft_bsl_stim(datapad,soi,'MISC004','MISC005',[-1 0],[0.25 1.25],'yes',[52 80],mGrad);
% 
% % induced
% [~, ~, bslind, rftind] = kd_fft_bsl_stim(datapad,soi,'MISC004','MISC005',[-1 0],[0.25 1.25],'no',[52 80],mGrad);
% 
% powchanev = rftev.powspctrm./bslev.powspctrm-1;
% powchanind = rftind.powspctrm./bslind.powspctrm-1;
% 
% fig = figure;
% subplot(221)
% plot(rftev.freq,rftev.powspctrm(1,:))
% title('evoked power spectrum')
% ylabel('power (T/m)²')
% subplot(222)
% plot(rftev.freq,powchanev(1,:))
% title('evoked power change')
% %ylabel('power change to baseline')
% subplot(223)
% plot(rftind.freq,rftind.powspctrm(1,:))
% title('induced power spectrum')
% ylabel('power (T/m)²')
% subplot(224)
% plot(rftind.freq,powchanind(1,:))
% title('induced power change')
% ylabel('power change to baseline')
% %legend(soicmb{:})
% 
% print(fig,fullfile(figpth,[subjcode,'_fft_',soicmb{1},'_sss']),'-dpng')
% %% Power spectral density
% 
% % identify reaction time
% rt = (cell2mat(cellfun(@length,data.time,'UniformOutput',false))-3.5*fs)./fs;
% 
% soi = {'MEG2032','MEG2033', 'MEG2042','MEG2043', 'MEG2112','MEG2113', 'MEG2122','MEG2123'};
% % select baseline
% cfg = [];
% cfg.channel = [soi, misc4,misc5];
% cfg.latency = [-1.25 -0.25];
% cfg.avgoverrpt = 'yes';
% BSLev = ft_selectdata(cfg,datapad);
% cfg.avgoverrpt = 'no';
% BSLind = ft_selectdata(cfg,datapad);
% 
% % select flicker evoked
% cfg.latency = [0.2 1.2];
% cfg.avgoverrpt = 'yes';
% FLICKev = ft_selectdata(cfg,datapad);
% cfg.avgoverrpt = 'no';
% FLICKind = ft_selectdata(cfg,datapad);
% 
% % trlflickind = {};
% % timeflickind = {};
% % trl = zeros(1,length(data.trial));
% % cfg.avgoverrpt = 'no';
% % 
% % % select flicker induced based on rt
% % for t = 1:length(rt)
% %     cfg.latency = [0.5 round(rt(t))];
% %     cfg.trials = trl;
% %     cfg.trials = t;
% %     FLICKind = ft_selectdata(cfg,data);
% %     trlflickind = [trlflickind, FLICKind.trial];
% %     timeflickind = [timeflickind,FLICKind.time];
% %     clear FLICKind
% % end
% % 
% % FLICKind = ft_selectdata(cfg,data);
% % FLICKind.time = timeflickind;
% % FLICKind.trial = trlflickind;
% 
% % pwelch PSD
% N = 512;
% for s = 1:length(soi)
%     [pbsl_ev(s,:), f] = pwelch(BSLev.trial{1}(s,:),hanning(N/2+1),N/2-10,[54:4:80],fs);
%     pfli_ev(s,:) = pwelch(FLICKev.trial{1}(s,:),hanning(N/2+1),N/2-10,[54:4:80],fs);
% end
% figure; plot(f,pfli_ev./pbsl_ev)
