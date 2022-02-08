%% VS + RFT
% PhD project 2

% Coherence per condition

function a_fun_coh_cond_new(s,condi,foi,fwdth,rej_sac)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','5 COH hilb','coh','sinusoid');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

d = dir(inpth);
d = {d.name};
folds = d(strncmp(d,'202',3));

if exist(fullfile(outpth,folds{s},strcat('coh_',strjoin(condi,'_'),'_freqw_',num2str(fwdth),'.mat')))
    error('coherence for subject and condition exists')
end

% list condition files that contain 'condi'
d = dir(fullfile(inpth,folds{s}));
if rej_sac
    d = dir(fullfile(inpth,folds{s},'rej_sac'));
    outpth = fullfile(outpth,'rej_sac');
end
d = {d.name};
%files = d(cell2mat(cellfun(@(x) ~isempty(x),regexp(d,condi),'UniformOutput',false)));
for c = 1:length(condi)
    cond_idx(c,:) = cell2mat(cellfun(@(x) ~isempty(x),regexp(d,condi(c)),'UniformOutput',false));
end
cond_idx = sum(cond_idx,1) == length(condi);
files = d(cond_idx);

load(fullfile(inpth,folds{s},files{1}));
data_load = data_trig;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,folds{s},files{f}));
    data_load = ft_appenddata([],data_load,data_trig);
end

data = data_load;

clear data_load

%% Coherence

% calculate coherence

% gradiometers
[coh.cohTgrad, powsd.psdmegTgrad, powsd.psdmiscTgrad, csd.csdTgrad] = kd_coh_hilb_fun(data,'diode T', 'MEGGRAD',foi, fwdth);
% coherence distractor
[coh.cohDgrad, powsd.psdmegDgrad, powsd.psdmiscDgrad, csd.csdDgrad] = kd_coh_hilb_fun(data,'diode D', 'MEGGRAD', foi, fwdth);

% number of trials
coh.ntrial = length(data.trial);


% magnetometers
[coh.cohTmag, powsd.psdmegTmag, powsd.psdmiscTmag, csd.csdTmag] = kd_coh_hilb_fun(data,'diode T', 'MEGMAG',foi, fwdth);
% coherence distractor
[coh.cohDmag, powsd.psdmegDmag, powsd.psdmiscDmag, csd.csdDmag] = kd_coh_hilb_fun(data,'diode D', 'MEGMAG', foi, fwdth);

mkdir(fullfile(outpth,folds{s}))
save(fullfile(outpth,folds{s},strcat('coh_',strjoin(condi,'_'),'_freqw_',num2str(fwdth),'.mat')),'coh','powsd','csd')




