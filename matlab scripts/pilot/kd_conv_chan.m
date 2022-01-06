%% Convolve hanning taper and data channel by channel

function compdat = kd_conv_chan(data,taper,N)

for c = 1:size(data,1)
    
    convdat = [];
    convdat = conv(data(c,:),taper);
    
    compdat(c,:) = convdat(ceil(N/2):length(convdat)-floor(N/2));
end