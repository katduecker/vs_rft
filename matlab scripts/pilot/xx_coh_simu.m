%% Coherence on simulated data

clear all; close all; clc; beep off;
mpth =  '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath(fullfile(mpth,'exgauss'))
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
ft_defaults;
mtlpth = fullfile(mpth,'matlab scripts');
f1 = 60;
f2 = 67;
fs = 1000;                          % sampling rate
bsl = 2.5;                          % baseline
numtrl = 100;                       % number of trials
numsens = 206;
idxdio1 = 205;
idxdio2 = 206;

load(fullfile(mtlpth,'templ_datastruct.mat'))% load template data structure

% adjust sample info
sampleinfo = repmat([1 8000],numtrl,1);
frsamp = [0:numtrl-1]'*8000;
sampleinfo = sampleinfo + repmat(frsamp,1,2);
datastruct.sampleinfo = sampleinfo;

rt = exgauss_rnd(1,0.25,1,numtrl);  % generate reaction time
rt(rt > 4) = 4;                     % cap at 4 seconds

% padded trials
trial = cell(1,numtrl);
time = cell(1,numtrl);
for t = 1:numtrl
    trial{t} = zeros(numsens,8*fs);
    time{t} = -2.5:1/fs:5.5-1/fs;
    timevec = 0:1/fs:rt(t);             % time vector current trial
    % diode signal
    trial{t}(idxdio1,1:bsl*fs+length(timevec)) = [zeros(1,bsl*fs), sin(2*pi*f1*timevec)];    % diode signal 1
    trial{t}(idxdio2,1:bsl*fs+length(timevec)) = [zeros(1,bsl*fs), sin(2*pi*f2*timevec)];    % diode signal 2

    % meg signal: both diodes + noise
    for c = 1:numsens-2
        noi = randn(1,bsl*fs+length(timevec));               % noise

        if c == find(strcmp(soi,datastruct.label))
            
            % RFT response in soi
            trial{t}(c,:) = trial{t}(c,:) + trial{t}(idxdio1,:) + trial{t}(idxdio2,:) + ...
                reshape([noi,zeros(1,8*fs-length(noi))],1,[]);
            
        else
            % rest noise
            trial{t}(c,:) = trial{t}(c,:) + reshape([noi,zeros(1,8*fs-length(noi))],1,[]);
        end
    end
end
datastruct.trial = trial;
datastruct.time = time;

% plot
timevec = -2.5:1/fs:5.5-1/fs;
subplot(411)
plot(timevec,squeeze(trial{1}(205,:)))
title('diode 1')
subplot(412)
plot(timevec,squeeze(trial{1}(206,:)))
title('diode 2')
subplot(413)
plot(timevec,squeeze(trial{1}(1,:)))
title('meg 1')
subplot(414)
plot(timevec,squeeze(trial{1}(find(strcmp(soi,datastruct.label)),:)))
title('meg soi')
xlabel('time (s)')
ylabel('frequency (Hz)')

% coherence fieldtrip
foi = 52:80;
frqwdth = 1;
[coh_spct, psd_meg, psd_misc,csd_meg_misc] = kd_coh_hilb(datastruct, foi, frqwdth, 'MISC004',8);


soiidx = find(strcmp(soi,datastruct.label));
fig = figure;
subplot(311)
imagesc(timevec,foi,squeeze(psd_misc(1,:,:)))
ylabel('frequency (Hz)')
title('PSD diode 1')
axis xy
subplot(312)
imagesc(timevec,foi,squeeze(psd_meg(1,:,:)))
ylabel('frequency (Hz)')
title('PSD meg random sensor')
axis xy
subplot(313)
imagesc(timevec,foi,squeeze(psd_meg(soiidx,:,:)))
ylabel('frequency (Hz)')
title('PSD meg soi')
axis xy
print(fig, fullfile(mpth,'pilot',['ft_simu_psd_frqwdth_',num2str(frqwdth)]),'-dpng')

fig = figure;
subplot(311)
imagesc(timevec,foi,squeeze(csd_meg_misc(205,:,:)))
ylabel('frequency (Hz)')
title('CSD diode 1')
axis xy
subplot(312)
imagesc(timevec,foi,squeeze(csd_meg_misc(1,:,:)))
ylabel('frequency (Hz)')
title('CSD meg random sensor')
axis xy
subplot(313)
imagesc(timevec,foi,squeeze(csd_meg_misc(soiidx,:,:)))
ylabel('frequency (Hz)')
title('CSD meg soi')
axis xy
print(fig, fullfile(mpth,'pilot',['ft_simu_csd_frqwdth_',num2str(frqwdth)]),'-dpng')

fig = figure;
subplot(211)
imagesc(timevec,foi,squeeze(coh_spct(1,:,:)))
ylabel('frequency (Hz)')
title('Coherence meg diode random sensor')
axis xy
subplot(212)
imagesc(timevec,foi,squeeze(coh_spct(soiidx,:,:)))
ylabel('frequency (Hz)')
title('Coherence meg diode soi')
axis xy
print(fig, fullfile(mpth,'pilot',['ft_simu_coh_frqwdth_',num2str(frqwdth)]),'-dpng')
