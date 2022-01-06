%% Pilot analysis: VS, alpha and RFT

%% Settings
clear all; close all; clc; beep off;
mpth = '/rds/projects/j/jenseno-visual-search-rft';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
ft_defaults;
subjcode = 'b57a';
soi = {'MEG2032','MEG2033','MEG2112','MEG2113'};
soicmb = {'MEG2032+2033','MEG2112+2113','MEG2042+2043','MEG1922+1923','MEG2342+2343'};
pthout = fullfile(mpth,'pilot','results',subjcode);

figpth = fullfile(mpth,'pilot','results','plots',subjcode);
mkdir(figpth)
fs = 1000;
load(fullfile(rsppth,[subjcode,'.mat']))
% load trigger
load(fullfile(cdpth,'trigdef.mat'))

d = dir(fullfile(megpth,subjcode));
files = {d.name};
files = files(strncmp(files,subjcode,4));

% just check color for now
files(cell2mat(cellfun(@isempty,strfind(files,'col'),'UniformOutput',false))) = [];

%% Load in data
% separately: misc004 picks up 60 Hz signal and misc005 picks up 60 Hz
% signal
%event = [];
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


% read in trials
for fl = 1:length(files)
    cfg = [];
    cfg.dataset  = fullfile(megpth,subjcode,files{fl});
    % get trial specifics
    cfg.path = cdpth;
    [trlsmp_misc004{fl}, trltrg_misc004{fl}, trl_misc004{fl},trlsmp_misc005{fl}, trltrg_misc005{fl}, trl_misc005{fl}] ...
        = kd_trlfun_trl_phd(cfg);
    cfg.trl = trl_misc004{fl};
    cfg.detrend = 'yes';
    %cfg.demean = 'yes';
    %cfg.baselinewindow = [-1.25 -.1];
    %cfg.hpfilter = 'yes';
    %cfg.hpfreq   = 40;
    block_misc004{fl} = ft_preprocessing(cfg);
    cfg.trl = trl_misc005{fl};
    block_misc005{fl} = ft_preprocessing(cfg);

end

% 1 block empty?
block_misc004(cell2mat(cellfun(@(x) isempty(x.trial),block_misc004,'UniformOutput',false)))= [];
block_misc005(cell2mat(cellfun(@(x) isempty(x.trial),block_misc005,'UniformOutput',false)))= [];

% concatenate
cfg = [];
cfg.keepsampleinfo = 'no';
data.misc004 = ft_appenddata(cfg,block_misc004{:});
data.misc005 = ft_appenddata(cfg,block_misc005{:}); 

% separate responses accordingly
[rspns.misc004, rspns.misc005] = kd_rspns_phd(subj);


%% Sanity checks

% plot trigger channel
trgchn = find(strcmp(data.misc004.label,'STI101'));

fig = figure;
for trl = 1:size(data.misc004.trial,2)
    curtrig = data.misc004.trial{trl}(trgchn,:);
    curtrig(curtrig > 36) = 50;
    curtrig(curtrig < 0) = 0;
    plot(data.misc004.time{trl},curtrig)
    hold on
end
print(fig,fullfile(figpth,'trigcheck'),'-dpng')
close all
% plot diode
fig = figure;
misc = find(strcmp(data.misc004.label,'MISC004'));
for trl = 1:size(data.misc004.trial,2)
    plot(data.misc004.time{trl},data.misc004.trial{trl}(misc,:))
 %  xlim([0 0.5])
    hold on
end
print(fig,fullfile(figpth,'misc004'),'-dpng')
close all
fig = figure;
misc = find(strcmp(data.misc005.label,'MISC005'));
for trl = 1:size(data.misc005.trial,2)
    plot(data.misc005.time{trl},data.misc005.trial{trl}(misc,:))
    xlim([0 0.5])
    hold on
end
print(fig,fullfile(figpth,'misc005'),'-dpng')
close all

% cut out long trials
l_004 = cell2mat(cellfun(@length,data.misc004.time,'UniformOutput',false));
l_005 = cell2mat(cellfun(@length,data.misc005.time,'UniformOutput',false));

% select data trials >= 1s
l_004trl = l_004 >= 4500;
l_005trl = l_005 >= 4500;

% time-lock analysis
cfg = [];
cfg.trials = l_004trl;
cfg.channel = {'MEG','MISC004','MISC005'};
data.misc004l = ft_selectdata(cfg,data.misc004);
cfg.trials = l_005trl;
data.misc005l = ft_selectdata(cfg,data.misc005);

cfg = [];
cfg.latency = [-2.5 2];
data.misc004l = ft_selectdata(cfg,data.misc004l);
data.misc005l = ft_selectdata(cfg,data.misc005l);

%% semi-automatic artefact rejection

cfg = [];
cfg.channel = 'MEGGRAD';
data.misc004meg = ft_selectdata(cfg,data.misc004l);
data.misc005meg = ft_selectdata(cfg,data.misc005l);

cfg = [];
cfg.method = 'summary';
cfg.layout = 'neuromag306planar.lay';
data.misc004cl = ft_rejectvisual(cfg,data.misc004meg);
data.misc005cl = ft_rejectvisual(cfg,data.misc005meg);

% find trials that are to be rejected
[~,rej4] = intersect(data.misc004l.sampleinfo(:,1),data.misc004cl.cfg.artfctdef.summary.artifact(:,1));
[~,rej5] = intersect(data.misc005l.sampleinfo(:,1),data.misc005cl.cfg.artfctdef.summary.artifact(:,1));

clear data.misc004meg data.misc005meg data.misc004cl data.misc005cl

cfg = [];
keeptrl = [1:length(data.misc004l.trial)];
keeptrl(rej4) = [];
cfg.trials = keeptrl;
data.misc004l = ft_selectdata(cfg,data.misc004l);

keeptrl = [1:length(data.misc005l.trial)];
keeptrl(rej5) = [];
cfg.trials = keeptrl;
data.misc005l = ft_selectdata(cfg,data.misc005l);

% append all trials for ERF
data.all = ft_appenddata([],data.misc004l,data.misc005l);

%% ERF
cfg = [];
cfg.keeptrials = 'no';
%ERF004 = ft_timelockanalysis(cfg,data.misc004l);
%ERF005 = ft_timelockanalysis(cfg,data.misc005l);
ERF = ft_timelockanalysis(cfg,data.all);
ERF.grad = mGrad;
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 40;
%cfg.demean = 'yes';
%cfg.baselinewindow = [-1.5 0];
ERFhp = ft_preprocessing(cfg,ERF);

cfg = [];
cfg.layout = 'neuromag306planar.lay';
cfg.baseline = [-1.5 -0.1];
cfg.xlim = [-1 1];
%ft_multiplotER(cfg,ERFhp);

for c = 1:length(soi)
    plot(ERFhp.time,ERFhp.avg(strcmp(ERFhp.label,soi{c}),:))
    xlim([-1,1])
    title(soi{c})
    xlabel('time (s)')
    ylabel('amplitude fT')
    print(fullfile(figpth,[subjcode,'_',soi{c},'_ERF_bp']),'-dpng')
   % pause
    close all
end


% TFR
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
cfg.foi = 0:1:100;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi));
cfg.toi = -1:0.05:.5;
cfg.keeptrials = 'no';
TFR = ft_freqanalysis(cfg,ERFhp);
TFR.grad = mGrad;

% ERF PSD
cfg = [];
cfg.latency = [0 1];
cfg.channel = {'MEG2112','MEG2113'};
ERFstim = ft_selectdata(cfg,ERF);
cfg.latency = [-1.25 -0.25];
ERFbsl = ft_selectdata(cfg,ERF);

% FFT ERF
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'hanning';
cfg.foilim = [40 100];
cfg.tapsmofrq = 1;
FFTstim = ft_freqanalysis(cfg,ERFstim);
FFTbsl = ft_freqanalysis(cfg,ERFbsl);

% combine planar
cfg = [];
cfg.method = 'sum';
FFTstim = ft_combineplanar(cfg,FFTstim);
FFTbsl = ft_combineplanar(cfg,FFTbsl);

% relative change
rc = FFTstim.powspctrm./FFTbsl.powspctrm -1;

% plot
plot(FFTstim.freq,rc)
xlabel('frequency (Hz)')
xlim([40 80])
ylabel('power change compared to bsl')
title('MEG2112+2113')
print(fullfile(figpth,[subjcode,'FFT_ERF_MEG2112+2113']),'-dpng')



%% plot

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);
% TFR004 = ft_combineplanar(cfg,TFR004);
% TFR005 = ft_combineplanar(cfg,TFR005);


close all;
cfg = [];
cfg.colorbar = 'yes';
cfg.baseline = [-1.5 -0.1];
cfg.baselinetype = 'relchange';
%cfg.zlim = [-2.5 2.5];
cfg.ylim = [55 70];
cfg.xlim = [-1 0.5];
cfg.layout = 'neuromag306cmb.lay';
%ft_multiplotTFR(cfg,TFR)
% 

for c = 1:length(soicmb)
    cfg.channel = soicmb{c};
    %cfg.zlim = [0 4]
    ft_singleplotTFR(cfg,TFR)
     print(fullfile(figpth,[subjcode,'TFR_pow_',soicmb{c}]),'-dpng')
   close all
end
% SOI = {};
% SOI_64 = {};
% for c = 2:length(TFR.label)
%     cfg.channel = TFR.label{c};
%     ft_singleplotTFR(cfg,TFR)
%     s = input('SOI? 1/0');
%     if s
%         SOI = [SOI, cfg.channel];
%     end
%     
%     s = input('SOI 64? 1/0');
%     if s
%         SOI_64 = [SOI_64, cfg.channel];
%     end
%     close all
% end
% 
% cfg.channel = SOI;
% ft_multiplotTFR(cfg,TFR)
% print(fullfile(figpth,'SOI_TFR_ERF'),'-dpng')
% 
% SOI60 = {};
% SOI67 = {};
% SOI6067 = {};
% for s = 1:length(SOI)
%     cfg.channel = SOI{s};
%     ft_singleplotTFR(cfg,TFR)
%     w = input('SOI? 1/2/3');
%     if w == 1
%         SOI60 = [SOI60, SOI{s}];
%     elseif w == 2
%         SOI67 = [SOI67, SOI{s}];
%     elseif w == 3
%         SOI6067 = [SOI6067, SOI{s}];
%     end
%     clear w
%     close all
% end
% 
% SOI67 = [SOI67, SOI60{5}];
% SOI60(5) = [];
% SOI = [SOI60,SOI67];
% 
% for s = 1:length(SOI60)
%     cfg.channel = SOI60{s};
%     cfg.xlim = [-1 .5];
%     cfg.ylim = [59 70];
%     subplot(3,2,s)
%     ft_singleplotTFR(cfg,TFR)
%     hold on
%     line(TFR.time,repmat(60,1,length(TFR.time)),'Color','black')
%     hold on
%     line(TFR.time,repmat(67,1,length(TFR.time)),'Color','black')
%     xlim([-1 0.5])
% end
% 
% cfg.channel = SOI60;
% ft_multiplotTFR(cfg,TFR);
% 
% save(fullfile(pthout,'TFR_ERF.mat'),'ERF','TFR','SOI','SOI60','SOI67','SOI6067','SOI_64')
% 
% savefig(fullfile(figpth,[subjcode,'TFR_pow_all.fig']))
% 
% soitfr = {'MEG1632+1633','MEG1912+1931','MEG2042+2043','MEG2142+2143','MEG2532+2533','MEG0132+0133'};
% save(fullfile(mpth,'pilot','results',[subjcode,'_ERF.mat']),'ERF','TFR','ERF004','ERF005','TFR004','TFR005','soitfr')
% close all
for soi = 1:length(soitfr)
    cfg.channel = soitfr{soi};
    ft_singleplotTFR(cfg,TFR)
    print(fullfile(figpth,[subjcode,'TFR_pow_',soitfr{soi}]),'-dpng')
    close all
end
clear ERF* TFR*

%% coherence
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'powandcsd';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
cfg.channelcmb = {'MEGGRAD','MISC004'};
%cfg.pad = 4;
%cfg.tapsmofrq = 2;
cfg.foi = 0:1:100;
cfg.taper = 'hanning';
cfg.t_ftimwin = 20./cfg.foi;
cfg.toi = -2:0.05:1.5;
cfg.keeptrials = 'yes';
TFR004 = ft_freqanalysis(cfg,data.misc004l);
TFR005 = ft_freqanalysis(cfg,data.misc005l);

TFR004.grad = mGrad;
TFR005.grad = mGrad;

% coherence
cfg = [];
cfg.method = 'coh';
cfg.channelcmb = {'MEGGRAD','MISC004'};
COH004 = ft_connectivityanalysis(cfg,TFR004);
COH005 = ft_connectivityanalysis(cfg,TFR005);
save(fullfile(mpth,'pilot','results',[subjcode,'_COH_misc004.mat']),'COH004')

% plot
cfg = [];
cfg.parameter = 'cohspctrm';
cfg.channel = 'MEGGRAD';
cfg.refchannel = 'MISC004';
cfg.lyout = 'neuromag306planar.lay';
cfg.colorbar = 'yes';
%cfg.zlim = [0 0.25];
cfg.ylim = [55 70];
ft_multiplotTFR(cfg,COH004)

soicoh = {'MEG2133','MEG0213','MEG1932','MEG2242','MEG2033','MEG2023'};
close all
for soi = 1:length(soicoh)
    cfg.channel = soicoh{soi};
    ft_singleplotTFR(cfg,COH004)
    print(fullfile(figpth,[subjcode,'COH004_',soicoh{soi}]),'-dpng')
    close all
end
hold on 
line(COH005.time,repmat(67,1,length(COH005.time)),'Color','black')
% all trials
% fund current set size
% for t = 5:size(trigdef,1)
%     trig32(t) = strcmp(trigdef{t,2}(3:4),num2str(32));
%     trig16(t) = strcmp(trigdef{t,2}(3:4),num2str(16));
% end
% ntrl = 0;

% select set sizes
% for fl = 1:size(block,2)
%     [~,cfg.trials] = find(ismember(trltrg{fl},[trigdef{logical(trig32)}]));
%     block32{fl} = ft_selectdata(cfg,block_keep{fl});
%     ntrl = ntrl + length(cfg.trials);
%     [~,cfg.trials] = find(ismember(trltrg{fl},[trigdef{logical(trig16)}]));
%     block16{fl} = ft_selectdata(cfg,block_keep{fl});
%     ntrl = ntrl + length(cfg.trials);
% 
% end
% 
% cfg = [];
% cfg.keepsampleinfo = 'no';
% data32 = ft_appenddata(cfg,block32{:});
% data16 = ft_appenddata(cfg,block16{:});

% cut trials based on shortest RT
% l32 = cell2mat(cellfun(@length,data32.time,'UniformOutput',false));
% l16 = cell2mat(cellfun(@length,data16.time,'UniformOutput',false));

% select data trials >= 1s
ltrl32 = l32 >= 2250;
ltrl16 = l16 >= 2250;

% time-lock analysis
cfg = [];
cfg.trials = ltrl32;
cfg.channel = {'MEG','MISC004','MISC005'};
data32cut = ft_selectdata(cfg,data32);
cfg.trials = ltrl16;
data16cut = ft_selectdata(cfg,data16);

cfg = [];
cfg.latency = [-1.5 1];
data32cut = ft_selectdata(cfg,data32cut);
data16cut = ft_selectdata(cfg,data16cut);

cfg = [];
ERF32 = ft_timelockanalysis(cfg,data32cut);
ERF16 = ft_timelockanalysis(cfg,data16cut);

% TFR on ERF
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'powandcsd';
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
cfg.channelcmb = {'MEGGRAD','MISC004'; 'MEG','MISC005'};

%cfg.channelcmb = {'MEGGRAD','MISC'};
cfg.pad = 3;
cfg.tapsmofrq = 2;
cfg.foi = 0:2:100;
cfg.taper = 'hanning';
cfg.t_ftimwin = ones(size(cfg.foi)) .* 0.5;
cfg.toi = -1.25:0.05:0.75;

TFR32 = ft_freqanalysis(cfg,ERF32);
TFR16 = ft_freqanalysis(cfg,ERF16);

cfg = [];
cfg.method = 'coh';
%cfg.channelcmb = {'MEG' 'MISC'};
COH32rft = ft_connectivityanalysis(cfg,TFR32);
COH16rft = ft_connectivityanalysis(cfg,TFR16);

% FFT
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.pad = 6;
cfg.channel     = {'MISC004','MISC005','MEG'};
cfg.taper       = 'hanning';
cfg.foilim      = [0 100];
cfg.tapsmofrq   = 2;
FREQ32y60       = ft_freqanalysis(cfg,data_misc004);
FREQ32t60       = ft_freqanalysis(cfg,data_misc005);

FREQ32rft       = ft_freqanalysis(cfg,data32);
FREQ16rft       = ft_freqanalysis(cfg,data16);
FREQ32bsl       = ft_freqanalysis(cfg,data32bsl);
FREQ16bsl       = ft_freqanalysis(cfg,data16bsl);

% coherence
cfg.output = 'powandcsd';
CSD32rft = ft_freqanalysis(cfg,data32);
CSD16rft = ft_freqanalysis(cfg,data16);



cfg.method = 'plv';
PLV32rft = ft_connectivityanalysis(cfg,CSD32rft);
PLV16rft = ft_connectivityanalysis(cfg,CSD16rft);


FREQ32y60.grad = mGrad;
FREQ32t60.grad = mGrad;

save(fullfile(mpth,'pilot','results',[subjcode,'_diodecheck.mat']),'FREQ32y60','FREQ32t60')

PLV32rft.grad = mGrad;
PLV16rft.grad = mGrad;
save(fullfile(mpth,'pilot','results',[subjcode,'_plv_setsize.mat']), 'PLV32rft','PLV16rft')

FREQ32rft.grad = mGrad;
FREQ16rft.grad = mGrad;
FREQ32bsl.grad = mGrad;
FREQ16bsl.grad = mGrad;
COH32rft.grad = mGrad;
COH16rft.grad = mGrad;
ERF32.grad = mGrad;
ERF16.grad = mGrad;
TFR32.grad = mGrad;
TFR16.grad = mGrad;

save(fullfile(mpth,'pilot','results',[subjcode,'_fft_setsize_short.mat']),'FREQ32rft','FREQ16rft','FREQ32bsl','FREQ16bsl')


save(fullfile(mpth,'pilot','results',[subjcode,'_coh_setsize.mat']),'ERF32','ERF16','TFR32','TFR16','COH32rft','COH16rft')

% plot
load(fullfile(mpth,'pilot','results',[subjcode,'_coh_setsize.mat']))

cfg = [];
ft_multiplotER(cfg,ERF32)


cfg = [];
cfg.method = 'sum';
TFRcmb32 = ft_combineplanar(cfg,TFRgr32);

cfg = [];
cfg.baseline = [-1.5 0];
cfg.zlim = [-3 3];
cfg.colorbar = 'yes';
%cfg.layout = 'neuromag306planar.lay';
cfg.baselinetype = 'relchange';
ft_multiplotTFR(cfg,TFRcmb32)

cfg = [];
cfg.parameter = 'cohspctrm';
cfg.refchannel = 'MISC004';
ft_multiplotTFR(cfg,COH32rft)



cfg.channel = {'MEG2042+2043'};
ft_singleplotTFR(cfg,TFRcmb32)
hold on 
line(TFRcmb32.time,repmat(60,1,length(TFRcmb32.time)),'Color','black')
hold on
line(TFRcmb32.time,repmat(67,1,length(TFRcmb32.time)),'Color','black')


cfg.channel = {'MEG2042'};

cfg.parameter = 'cohspctrm';
cfg.refchannel = 'MISC004';
ft_singleplotTFR(cfg,COHgr32)

%% check color amplitude