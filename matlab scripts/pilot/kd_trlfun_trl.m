
function trl = kd_trlfun_trl_phd(cfg)

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);

% trigger + button presses
allval = [event(find(strcmp('STI101',{event.type}))).value];   % these are all values, extract button press
allsmp = [event(find(strcmp('STI101',{event.type}))).sample];

% trial structure: condition trigger - 2 - button
a = 1;
% trlsmp = [];
% trltrg = [];
trl    = [];
x = [];
load(fullfile(cfg.path,'trigdef.mat'))
while a < length(allval)
    if allval(a) >= 5 && allval(a) <= trigdef{end,1}
        bslsmp = allsmp(a);
        if allval(a+1) == 2
            sdsmp = allsmp(a+1);
        else
            sdsmp = NaN;
            x = [x,a];
        end
        if allval(a+2) > trigdef{end,1}
            btsmp = allsmp(a+2);
            
            % delete 6 -> just because trigger 2 wasn't resent
        elseif allval(a+2) == 4 || allval(a+2) == 6
            % if no button: nan
            btsmp = NaN;
            x = [x,a];
        else
            % if no button: nan
            btsmp = NaN;
            x = [x,a];
        end
        
        % shouldn't happen
        if isnan(sdsmp)
            sdsmp = bslsmp + hdr.Fs*1.5;
        end
        
        if isnan(btsmp)
            btsmp = sdsmp + hdr.Fs*4;
        end
        
        % allow time before baseline for filter ringing...
        if bslsmp-hdr.Fs > 0
%             trlsmp = [trlsmp;bslsmp sdsmp btsmp];
%             trltrg = [trltrg, allval(a)];
            trl    = [trl; bslsmp-hdr.Fs btsmp+hdr.Fs -hdr.Fs*1.5];
        end
        a = a + 3;
        clear bslsmp sdsmp btsmp
    else 
        a = a + 1;
    end
    
end

