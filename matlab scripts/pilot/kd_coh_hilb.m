%% Coherence with hilbert transform

% Inputs: 
% - data: MEG data containing SOI and diodes
% - foi: frequencies of interest
% - frqwdth: width of passband
% - diodelabel: label of photodiode

% Output:
% coh_spect: coherence spectrum

function [coh_spct, psd_meg, psd_misc,csd_meg_misc] = kd_coh_hilb(data, diodedata,foi, frqwdth)

cfg = [];
cfg.bpfilter    = 'yes';
cfg.hilbert     = 'complex';
cfg.keeptrials  = 'yes';
%cfg.bpfilttype  = 'fir';
%cfg.plotfiltresp = 'yes';
%cfg.bpfiltdir = 'onepass-zerophase';
%cfg.trials = 1;
%cfg.padding = padl;
%cfg.bpfiltord = N;
%idxmisc = find(strcmp(data.label,diodelabel));
%padl = 8;       % length of padding in sec

% loop over frequencies and do hilbert transform of data
for f = 1:length(foi)
    % bp filter
    cfg.bpfreq = [foi(f)-frqwdth foi(f)+frqwdth];
    data_foi = ft_preprocessing(cfg,data);
    diode_foi = ft_preprocessing(cfg,diodedata);
    
    % init matrices with zero padding
    meg_pad = zeros(length(data.trial),length(data.label),length(data.time{1}));
    misc_pad = meg_pad;
    csd_pad = meg_pad;
    
    % for each channel
         
        % loop over trials
        for t = 1:length(data.trial)
            % magnitude per trial
            
            % meg in soi
            meg  = data_foi.trial{t};
            % diode
            misc = repmat(diode_foi.trial{t},size(meg,1),1);
            
            % power spectral density in this trial
            psd_meg_trl = meg.*conj(meg);
            psd_misc_trl = misc.*conj(misc);
            
            % cross-spectral density in this trial
            csd_meg_misc_trl = meg .* conj(misc);
            
            % fill into trialxchannelxtime matrix
            meg_pad(t,:,:) = psd_meg_trl;
            misc_pad(t,:,:) = psd_misc_trl;
            csd_pad(t,:,:) = csd_meg_misc_trl;
            
        end
    
    mean_meg_pad = squeeze(nanmean(meg_pad,1));
    mean_misc_pad = squeeze(nanmean(misc_pad,1));
    mean_csd_pad = squeeze(nanmean(csd_pad,1));
    
    % match fieldtrip format for plotting
    for c = 1:size(meg_pad,2)
        % psd: average over trials
        psd_meg(c,f,:) = mean_meg_pad(c,:);
        % average psd over trials -> copy such that size of psd of misc is the same
        % as psd of meg
        psd_misc(c,f,:) = mean_misc_pad(c,:);
        
        %csd by hand
        csd_meg_misc(c,f,:) = abs(mean_csd_pad(c,:)).^2;
    end
        
end
        
coh_spct = csd_meg_misc./(psd_meg.*psd_misc);           % coherence spectrum

% % only select -1.5 2 seconds (this is a bit of a stupid way of doing it...)
% [~,ilow] = min(abs(data.time{1} + 1.5));
% [~,iup] = min(abs(data.time{1} - 2));
% 
% psd_meg = psd_meg(:,:,ilow:iup);
% psd_misc = psd_misc(:,:,ilow:iup);
% csd_meg_misc = csd_meg_misc(:,:,ilow:iup);
% coh_spct = coh_spct(:,:,ilow:iup);
