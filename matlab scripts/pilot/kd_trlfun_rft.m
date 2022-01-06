function trl= kd_trlfun_trl(cfg)

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);

% trigger events
trgval = [event(find(strcmp('Trigger',{event.type}))).value];   % values
trgsmp = [event(find(strcmp('Trigger',{event.type}))).sample];  % samples

% button presses
allval = [event(find(strcmp('STI101',{event.type}))).value];   % these are all values, extract button press
allsmp = [event(find(strcmp('STI101',{event.type}))).sample];
% find first trigger after value 2 (trial beginning)
trlbidx = find(allval == 2);

btnidx = find(allval > 5000);
prsidx = [];
for t = 1:length(trlbidx)
    p = btnidx(btnidx - trlbidx(t) == 1);
    prsidx = [prsidx,p];
end

btnval = allval(prsidx);
btnsmp = allsmp(prsidx)';
%trlval = allval(prsidx-1) % sanity check
trlsmp = allsmp(prsidx-1)';

trl = [trlsmp btnsmp+hdr.Fs/2 zeros(size(trlsmp))];
rt = btnsmp' -trlsmp';
