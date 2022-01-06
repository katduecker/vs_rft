% Coherence function on simulated data

function [coh_spct, psd_meg, psd_misc,csd_meg_misc] = kd_coh_hilb_mtlb_simu(data, foi, freqwidth, N,)

idxmisc = strcmp(data.label,diodelabel);          % diode index

% loop over frequencies
for f = 1:length(foi)
    disp(['frequency ', num2str(foi(f))])

    % init matrices with zero padding
    meg_pad = zeros(length(data.trial),length(data.label),padl*data.fsample);
    misc_pad = meg_pad;
    csd_pad = meg_pad;
    
    % loop over trials
    for t = 1:length(data.trial)
        
        meg = [];
        misc = [];
        
        filtdat = kd_bp_filt('bp',data.trial{t},foi(f),N,freqwidth,fs);
        
        meg = hilbert(filtdat);
        misc = repmat(hilbert(filtdat(idxmisc,:)),size(meg,1),1);
        
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
        
coh_spct = csd_meg_misc./(psd_meg.*psd_misc);