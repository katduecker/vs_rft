%% Function: BP filter using HP and LP filter

% Inputs
% - data: to be filtered data
% - foi: frequency of interest single frequency
% - N: filter order/length of filter
% - frqwdth: transition band
% - fs: sampling rate

% Output:
% filtdat: bandpass filtered data

function [filtdat] = kd_bp_filt(fltp,data,foi,N,frqwdth,fs)

if strcmp(fltp,'bp')
    dbp = designfilt('bandpassiir','FilterOrder',N,'HalfPowerFrequency1',...
        foi-frqwdth,'HalfPowerFrequency2',foi+frqwdth,'SampleRate',1000);
    %[b,a] = butter(N/2,[foi-frqwdth foi+frqwdth]./(fs/2));
    for c = 1:size(data,1)
        filtdat(c,:) = filtfilt(dbp,data(c,:));
    end
elseif strcmp(fltp,'hplp')
    dh = designfilt('highpassfir','FilterOrder',N,'PassbandFrequency',...
        foi,'StopbandFrequency',foi-frqwdth,'SampleRate',fs);
    dl = designfilt('lowpassfir','FilterOrder',N,'PassbandFrequency',...
        foi,'StopbandFrequency',foi+frqwdth,'SampleRate',fs);
    %fvtool(dl)
    
    filtdat = filtfilt(dh,data');
    filtdat = filtfilt(dl,filtdat)';
    %     for c = 1:size(data,1)
    %         filtdat(c,:) = filtfilt(dh,data(c,:));
    %         filtdat(c,:) = filtfilt(dl,filtdat(c,:));
    %     end
end