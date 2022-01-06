
function d2_coh_contrasts(s,groupcode)
% s:         subject ID
% groupcode: trigger code common in trials 

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');

ldiode = {'MISC004','MISC005'};
% matlab scripts path
mtpth = fullfile(pth,'matlab scripts','tfrs');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb');
cohsoipth = fullfile(cohpth,'soi');
cohdatpth = fullfile(cohpth,'data sep');
cohrespth = fullfile(cohpth,'coh');
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft','fieldtrip'))
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

d = dir(fullfile(cohdatpth,subjfolds{s}));
folds = {d.name};
datfolds = folds(strncmp(folds,'data',4));

% find files that belong to current group
idx = cellfun(@(x) regexp(x,'6067'),datfolds,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
datfolds = datfolds(idxx);

load(fullfile(cohdatpth,subjfolds{s},datfolds{1}));

data_load = data;

for d = 2:length(datfolds)
    load(fullfile(cohdatpth,subjfolds{s},datfolds{d}));
    data_load = ft_appenddata([],data_load,data);
    clear data
end

data = data_load;
clear data_load
% % check resampling effect
% cfg = [];
% cfg.method ='mtmfft';
% cfg.taper = 'hanning';
% cfg.output = 'pow';
% cfg.channel = 'MEGGRAD';
% cfg.pad = 'nextpow2';
% freqdat = ft_freqanalysis(cfg,data);
% 
% plot(freqdat.freq,freqdat.powspctrm)

%% Coherence
foi = 52:76;
frqwdth = 2;

% Target freq 60 Hz
% gradiometers
[coh.cohTgrad, powsd.psdmegTgrad, powsd.psdmiscTgrad, csd.csdTgrad] = kd_coh_hilb_fun(data,'diode T', 'MEGGRAD',foi, frqwdth);
% coherence distractor
[coh.cohDgrad, powsd.psdmegDgrad, powsd.psdmiscDgrad, csd.csdDgrad] = kd_coh_hilb_fun(data,'diode D', 'MEGGRAD', foi, frqwdth);

% magnetometers
[coh.cohTmag, powsd.psdmegTmag, powsd.psdmiscTmag, csd.csdTmag] = kd_coh_hilb_fun(data,'diode T', 'MEGMAG',foi, frqwdth);
% coherence distractor
[coh.cohDmag, powsd.psdmegDmag, powsd.psdmiscDmag, csd.csdDmag] = kd_coh_hilb_fun(data,'diode D', 'MEGMAG', foi, frqwdth);


mkdir(fullfile(cohrespth,subjfolds{s}))
save(fullfile(cohrespth,subjfolds{s},['coh_',groupcode,'_freqw_',num2str(frqwdth),'.mat']),'coh','powsd','csd')
