%% VS + RFT
% PhD project 2

% Coherence per condition

function a_fun_coh_cond(s,condi,foi,fwdth)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions');
outpth = fullfile(pth,'results','meg','5 COH hilb','coh');
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
d = dir(inpth);
d = {d.name};
folds = d(strncmp(d,'202',3));

% list condition files that contain 'condi' (input)
d = dir(fullfile(inpth,folds{s}));
d = {d.name};
files = d(cell2mat(cellfun(@(x) ~isempty(x),regexp(d,condi),'UniformOutput',false)));

load(fullfile(inpth,folds{s},files{1}));
data_load = data;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,folds{s},files{f}));
    data_load = ft_appenddata([],data_load,data);
end

data = data_load;

clear data_load

%% Coherence

% calculate coherence

% gradiometers
[coh.cohTgrad, powsd.psdmegTgrad, powsd.psdmiscTgrad, csd.csdTgrad] = kd_coh_hilb_fun(data,'diode T', 'MEGGRAD',foi, fwdth);
% coherence distractor
[coh.cohDgrad, powsd.psdmegDgrad, powsd.psdmiscDgrad, csd.csdDgrad] = kd_coh_hilb_fun(data,'diode D', 'MEGGRAD', foi, fwdth);

% magnetometers
[coh.cohTmag, powsd.psdmegTmag, powsd.psdmiscTmag, csd.csdTmag] = kd_coh_hilb_fun(data,'diode T', 'MEGMAG',foi, fwdth);
% coherence distractor
[coh.cohDmag, powsd.psdmegDmag, powsd.psdmiscDmag, csd.csdDmag] = kd_coh_hilb_fun(data,'diode D', 'MEGMAG', foi, fwdth);

mkdir(fullfile(outpth,folds{s}))
save(fullfile(outpth,folds{s},['coh_',condi,'_freqw_',num2str(fwdth),'.mat']),'coh','powsd','csd')




