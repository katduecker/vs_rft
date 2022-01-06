function trl = kd_trlfun(cfg)

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
% clean trlidx
trlidx_cl = [];
for t = 1:length(trlbidx)
    p = btnidx(btnidx - trlbidx(t) == 1);
    if p
        prsidx = [prsidx,p];
        trlidx_cl = [trlidx_cl,trlbidx(t)];
    end
end

btnval = allval(prsidx);
btnsmp = allsmp(prsidx)';
%trlval = allval(prsidx-1) % sanity check
trlsmp = allsmp(trlidx_cl)';

trl = [trlsmp-hdr.Fs*2.5 btnsmp+hdr.Fs -ones(size(trlsmp)).*(hdr.Fs*2.5)];
if find(trl(:,1) < 0)
    trl(find(trl(:,1) < 0),:) = [];
end
rt = btnsmp' -trlsmp';


