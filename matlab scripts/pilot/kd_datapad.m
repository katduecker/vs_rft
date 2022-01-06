
function datapad = kd_datapad(data,fs,padl,time)

% pad data to 8 sec
datapad = data;
% change time vector
datapad.time(:) = {time};
padlength = padl*fs;
% apply zero padding
% target tagged at 60
for t = 1:length(datapad.trial)
    pad = zeros(length(data.label),padlength-size(data.trial{t},2));
    datapad.trial{t} = horzcat(datapad.trial{t},pad);
end