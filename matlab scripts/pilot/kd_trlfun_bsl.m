function trl= kd_trlfun_bsl(cfg)

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);

% trigger events
trgval = [event(find(strcmp('Trigger',{event.type}))).value];   % values
trgsmp = [event(find(strcmp('Trigger',{event.type}))).sample];  % samples

% find trigger of trial beginning: 2
trlbidx = find(trgval == 2);
% baseline
bslidx = find(((trgval > 4) + (trgval<5000)) == 2);

bsllength = trgsmp(trlbidx) - trgsmp(bslidx);              % check if this is as long as bsl! -> yes!!
trlsmp = trgsmp(trlbidx)';

trl = [trlsmp-hdr.Fs*1.5 trlsmp -ones(size(trlsmp)).*(hdr.Fs*1.5)];

