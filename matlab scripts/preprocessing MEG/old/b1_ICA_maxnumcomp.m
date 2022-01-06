%% VS + RFT
% PhD project 2

% artefact suppression using ICA
% step 1: find number of principal components -> maximum number of ICA
% comps
% [c] Katharina Duecker

clear all; close all; clc; beep off
% define paths
pth = 'W:\Visual Search RFT';
addpath(fullfile('matlab scripts','edf-converter-master')); %edf2mat converter
dtpth = fullfile(pth,'results','meg', '1 maxfilter');
addpath(fullfile('W:\','fieldtrip'))            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

%numcomp = zeros(1,length(subjfolds));            % store maximum number of components for each subject
load('ICA_num_comp.mat')
%% Maximum number of components per subject
for s = 20%length(subjfolds)
    % load in trial structure
    load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'))
    
    % trial structure to load in trl
    for p = 1%:2%length(meginfo.alltrl_bl)
        trlstruct{p} = [meginfo.alltrl_bl{p}(:,2)-fs,meginfo.alltrl_bl{p}(:,4)+fs,zeros(length(meginfo.alltrl_bl{p}),1)-1.5];
        trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;

    end
    
    % load in maxfiltered data using trial structure
    
    % list fif files
    d = dir(fullfile(dtpth,subjfolds{s}));
    f = {d.name};
    % find fif file
    idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    f = f(idxx);
    
    cfg = [];
    for p = 1%:2%length(f)
        cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
        cfg.detrend = 'yes';
        cfg.trl = trlstruct{p};
        cfg.channel = 'MEG';
        % load in data for this part
        dtprt{p} = ft_preprocessing(cfg);
    end
    
    data = ft_appenddata([],dtprt{:});
    
    % fix sample info for downsampling  
    sinfo = zeros(length(data.trial),1);
    % first trial samples: 0 to length of data in samples
    sinfo(1,2) = size(data.trial{1},2);
    for t = 2:length(data.trial)
        % next trial: one sample after end of previous
        sinfo(t,1) = sinfo(t-1,2) + 1;
        % start sample + length of trial - 1
        sinfo(t,2) = sinfo(t,1) + size(data.trial{t},2);
    end
    
    data.sampleinfo = sinfo;
    
    % resample to 250 Hz (1/4 of sampling rate)
%     
%     cfg = [];
%     cfg.resamplefs = 250;
%     dataresamp = ft_resampledata(cfg,data);
    
    datapad = data;
    clear data
    % check rank of the data
    % zero-pad data
    padlength = 10*fs;
    
    % change sampleinfo
    datapad.sampleinfo = horzcat([1:padlength:padlength*length(datapad.trial)]',[1:padlength:padlength*length(datapad.trial)]'+padlength-1);
    
    timevec = linspace(0,padlength/fs,padlength);
    % apply zero padding
    for t = 1:length(datapad.trial)
        pad = zeros(length(datapad.label),padlength-size(datapad.trial{t},2));
        datapad.trial{t} = horzcat(datapad.trial{t},pad);
        datapad.time{t} = timevec;
    end
    
    % covariance matrix
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'all'; %can be all, or time window of interest
    cfg.keeptrials = 'no';

    data_cov = ft_timelockanalysis(cfg, datapad);
    clear datapad
    % svd
    [~,svd_decomp,~] = svd(data_cov.cov);
    
    % plot
    figure;
    semilogy(diag(svd_decomp),'o-');
    % find cut off based on 10 smallest gradients
    x = gradient(log(diag(svd_decomp)));
    % difference of gradients -> find sudden drop
    y = diff(x);
    % difference minimum >2*std
%     localminy = find(islocalmin(y,'MinProminence',2*std(y)));
%     cutoff = localminy(1);
%     prompt = input(['max number of comps ? y'], 's');
%     
%     if strcmp(prompt,'y')
%         
%         numcomp(s) = cutoff;
%     else
%         numcomp(s) = input('cutoff?');
%     end
    numcomp(s) = input(['max number of comps ?']);
    close all
    clear data* x y cutoff localminy prompt
end
 
save('ICA_num_comp.mat','numcomp')