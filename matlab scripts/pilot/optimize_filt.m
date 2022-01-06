%% Optimize filter

N = 4;
foi = 50;%:1:80;
frqwdth = 2;
fs = 1000;
for f = 1:length(foi)
    [A,B,C,D] = butter(N/2,[foi-frqwdth foi+frqwdth]./(1000/2));
    sos = ss2sos(A,B,C,D);
    dbp = designfilt('bandpassiir','FilterOrder',N,'HalfPowerFrequency1',...
        foi-frqwdth,'HalfPowerFrequency2',foi+frqwdth,'SampleRate',1000);
    %fvt = fvtool(sos,dbp,'Fs', 1000);
    %legend(fvt,'butter','designfilt')
    
    figure;
    t = linspace(-2.5,5.5,8000);
    sig = [zeros(1,length(t)/2).*0.5,sin(2*pi*50*t(length(t)/2 + 1:end))];
    sig = repmat(sig,206,1);
    subplot(311)
    plot(t,sig)
    %xlim([0 0.1])
    title('original signal')
    sig = sig + randn(size(sig));
    subplot(312)
    plot(t,sig)
    %xlim([0 0.1])
    title('signal + noise')
    subplot(313)
    filtsig = kd_bp_filt('bp',sig,foi,N,frqwdth,fs);
    plot(t,filtsig)
    %xlim([0 0.1])
    title('filtered signal')
%     dh = designfilt('highpassfir','FilterOrder',N,'PassbandFrequency',...
%         foi(f),'StopbandFrequency',foi(f)-frqwdth,'SampleRate',fs);
%     fvtool(dh)
%     dl = designfilt('lowpassfir','FilterOrder',N,'PassbandFrequency',...
%         foi(f),'StopbandFrequency',foi(f)+frqwdth,'SampleRate',fs);
%     fvtool(dl)
    
    % impulse response
    imp = zeros(1,8001);
    imp(4000) = 1;
    impresp = kd_bp_filt('bp',imp,foi,N,frqwdth,fs);
  
    fig = figure;
    subplot(211)
    plot(impresp)
    xlabel('time (s)')
    title('impulse response IIR filter')
    xlim([3000 5000])
    % plot frequency response
    freqimp = fft(impresp);
    P2 = abs(freqimp/length(impresp));
    P1 = P2(1:length(impresp)/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    freq = fs*(0:(length(impresp)/2))/length(impresp);
    subplot(212)
    plot(freq,P1)
    hold on
    line(freq,repmat(max(P1)/2,1,length(freq)),'Color','k','LineStyle','-.')
    xlabel('frequency (Hz)')
    title('frequency response')
    xlim([0 100])
    pause
    close all
end
print(fig, fullfile('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/pilot',['FIR_order_',num2str(N),'_49_50Hz']),'-dpng')