%% VS + RFT
% PhD project 2

% replace diode signal with sinusoid
% ranging from 0 to 1, amplitude 0.5, mean 0.5

% [c] Katharina Duecker

function data = kd_replace_diode_sinu(data,diode_idx,rft_freq,phshft,start_bsl,fs)

% inputs
% - data: data structure (fieldtrip)
% - diode_idx: diode indices
% - f: frequencies
% - phshft: phase shift of sinusoid
% - start_bsl: time point of bsl start relative to display onset in s

for t = 1:length(data.trial)
    for di = 1:length(diode_idx)
        grad_diode = round(gradient(data.trial{t}(diode_idx(di),abs(start_bsl*fs):end)),3);
        diode_delay = find(grad_diode>0,1);
        if phshft(di) > 0
            diode_delay = find(grad_diode>0,1)+1;
        end
        tvec = 0:1/fs:data.time{1}(end)-diode_delay/fs+1/fs;
        rft_sig = [zeros(1,abs(start_bsl*fs)+diode_delay-1),sin(2*pi*rft_freq(di)*tvec+phshft(di)).*0.5+0.5];
        
        data.trial{t}(diode_idx(di),:) = rft_sig;
    end
end