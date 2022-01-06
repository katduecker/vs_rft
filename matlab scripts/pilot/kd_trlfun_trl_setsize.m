function trl = kd_trlfun_trl(cfg)

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);
load(fullfile(cfg.path,'trigdef.mat'))
% set size of interest trigger
ssoi = zeros(size(trigdef,2),1);
% find current set size
for t = 5:size(trigdef,1)
    ssoi(t) = strcmp(trigdef{t,2}(3:4),num2str(cfg.setsize));
end

% trigger + button presses
allval = [event(find(strcmp('STI101',{event.type}))).value];   % these are all values, extract button press
allsmp = [event(find(strcmp('STI101',{event.type}))).sample];

% trial structure: condition trigger - 2 - button
a = 1;
trlsmp = [];
trltrg = [];
trl    = [];
while a < length(allval)
    if find(allval(a) == [trigdef{logical(ssoi),1}])
        bslsmp = allsmp(a);
        if allval(a+1) == 2
            sdsmp = allsmp(a+1);
        else
            sdsmp = NaN;
        end
        if allval(a+2) > trigdef{end,1}
            btsmp = allsmp(a+2);
        else
            % if no button: nan
            btsmp = NaN;
        end
        trlsmp = [trlsmp;bslsmp sdsmp btsmp];
        trltrg = [trltrg, allval(a)];
        if isnan(btsmp)
            btsmp = sdsmp + hdr.Fs*4;
        end
        % shouldn't happen
        if isnan(sdsmp)
            sdsmp = bslsmp + hdr.Fs*1.5;
        end
        
        btsmp = btsmp + hdr.Fs;

        trl    = [trl; sdsmp-hdr.Fs*2.5 btsmp -hdr.Fs*2.5];
        a = a + 3;
        clear bslsmp sdsmp btsmp
    else 
        a = a + 1;
    end
    
end

