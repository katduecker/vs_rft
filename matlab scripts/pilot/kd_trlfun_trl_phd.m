% Function separates data into 2 data sets
% data set 1: MISC004 picks up 60 Hz RFT, MISC005 67 Hz
% 2: MISC004 67 Hz, MSC005 60 Hz

function [trlsmp_misc004, trltrg_misc004, trl_misc004,trlsmp_misc005,trltrg_misc005,trl_misc005] = kd_trlfun_trl_phd(cfg)

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);

% trigger + button presses
allval = [event(find(strcmp('STI101',{event.type}))).value];   % these are all values, extract button press
allsmp = [event(find(strcmp('STI101',{event.type}))).sample];

% trial structure: condition trigger - 2 - button
a = 1;
trlsmp_misc004 = [];
trltrg_misc004 = [];
trl_misc004    = [];

trlsmp_misc005 = [];
trltrg_misc005 = [];
trl_misc005    = [];

x = [];
load(fullfile(cfg.path,'trigdef.mat'))

% trigger of interest

% color 1 is target & 60 Hz
f60_misc004_t = 4 + find((cell2mat(cellfun(@(x) strcmp(x(7),'1'),trigdef(5:end,2),'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(9:10),'60'),trigdef(5:end,2),'UniformOutput',false))) == 2);
% color 1 is distractor & 60 Hz
f60_misc004_d = 4 + find((cell2mat(cellfun(@(x) strcmp(x(8),'1'),trigdef(5:end,2),'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(11:12),'60'),trigdef(5:end,2),'UniformOutput',false))) == 2);
% trigger values
trig60_misc004 = [trigdef{[f60_misc004_t;f60_misc004_d],1}];

% color 2 is target & 60 Hz
f60_misc005_t = 4 + find((cell2mat(cellfun(@(x) strcmp(x(7),'2'),trigdef(5:end,2),'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(9:10),'60'),trigdef(5:end,2),'UniformOutput',false))) == 2);
% color 2 is distractor & 60 Hz
f60_misc005_d = 4 + find((cell2mat(cellfun(@(x) strcmp(x(8),'2'),trigdef(5:end,2),'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(11:12),'60'),trigdef(5:end,2),'UniformOutput',false))) == 2);
trig60_misc005 = [trigdef{[f60_misc005_t;f60_misc005_d],1}];

clear f60*

t = 0;
while a+2 < length(allval)
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
        
        % add one second of data before and after trial
        bslsmp = bslsmp - hdr.Fs;
        btsmp = btsmp + hdr.Fs;
        % allow time before baseline for filter ringing...
        if bslsmp-hdr.Fs > 0
            if find(trig60_misc004 == allval(a))
                trlsmp_misc004 = [trlsmp_misc004;bslsmp sdsmp btsmp];
                trltrg_misc004 = [trltrg_misc004, allval(a)];
                trl_misc004    = [trl_misc004; bslsmp btsmp -hdr.Fs*2.5];
            elseif find(trig60_misc005 == allval(a))
                trlsmp_misc005 = [trlsmp_misc005;bslsmp sdsmp btsmp];
                trltrg_misc005 = [trltrg_misc005, allval(a)];
                trl_misc005    = [trl_misc005; bslsmp btsmp -hdr.Fs*2.5];
            end
        end
        a = a + 3;
        clear bslsmp sdsmp btsmp
    else 
        a = a + 1;
    end
    
end
