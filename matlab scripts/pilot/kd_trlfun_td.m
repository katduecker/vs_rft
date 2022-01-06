%% Trigger based on Target

function [trl60_misc4, trl60_misc5, trl67_misc4, trl67_misc5] = kd_trlfun_td(cfg)

hdr     = ft_read_header(cfg.dataset);
event   = ft_read_event(cfg.dataset);

% trigger + button presses
allval = [event(find(strcmp('STI101',{event.type}))).value];   % these are all values, extract button press
allsmp = [event(find(strcmp('STI101',{event.type}))).sample];

load(fullfile(cfg.path,'trigdef.mat'))

% target tagged
if strcmp(cfg.tag, 't')
    % trigger of interest
    % color 1 is target & 60 Hz
    f60_misc004 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(7),'1'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(9:10),'60'),trigdef(5:end,2),...
        'UniformOutput',false))) == 2);
    % color 2 is target & 60 Hz
    f60_misc005 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(7),'2'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(9:10),'60'),trigdef(5:end,2),...
        'UniformOutput',false))) == 2);
    
    trig60_misc4 = [trigdef{f60_misc004,1}];
    trig60_misc5 = [trigdef{f60_misc005,1}];
    
    % color 1 is target & 67 Hz
    f67_misc004 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(7),'1'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(9:10),'67'),...
        trigdef(5:end,2),'UniformOutput',false))) == 2);
    % color 2 is target & 67 Hz
    f67_misc005 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(7),'2'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(9:10),'67'),trigdef(5:end,2),...
        'UniformOutput',false))) == 2);
    
    trig67_misc4 = [trigdef{f67_misc004,1}];
    trig67_misc5 = [trigdef{f67_misc005,1}];
elseif strcmp(cfg.tag, 'd')
    % color 1 is distractor & 60 Hz
    f60_misc004 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(8),'1'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(11:12),'60'),trigdef(5:end,2),...
        'UniformOutput',false))) == 2);
    % color 2 is distractor & 60 Hz
    f60_misc005 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(8),'2'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(11:12),'60'),trigdef(5:end,2),...
        'UniformOutput',false))) == 2);
    
    trig60_misc4 = [trigdef{f60_misc004,1}];
    trig60_misc5 = [trigdef{f60_misc005,1}];
    
    % color 1 is distractor & 67 Hz
    f67_misc004 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(8),'1'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(11:12),'67'),...
        trigdef(5:end,2),'UniformOutput',false))) == 2);
    % color 2 is distractor & 67 Hz
    f67_misc005 = 4 + find((cell2mat(cellfun(@(x) strcmp(x(8),'2'),trigdef(5:end,2),...
        'UniformOutput',false)) + cell2mat(cellfun(@(x) strcmp(x(11:12),'67'),trigdef(5:end,2),...
        'UniformOutput',false))) == 2);
    
    trig67_misc4 = [trigdef{f67_misc004,1}];
    trig67_misc5 = [trigdef{f67_misc005,1}];
end

trl60_misc4 = [];
trl60_misc5 = [];
trl67_misc4 = [];
trl67_misc5 = [];

% define trials
a = 1;          % while loop counter
x = [];         % safe faulty trials
% trial structure: condition trigger - 2 - button

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
        btsmp = btsmp + 2*hdr.Fs;
        
        % allow time before baseline for filter ringing...
        if bslsmp-hdr.Fs > 0
            % misc4 tagged at 60
            if find(trig60_misc4 == allval(a))
                trl60_misc4 = [trl60_misc4;bslsmp btsmp -hdr.Fs*2.5];
            
            % if misc5 tagged at 60
            elseif find(trig60_misc5 == allval(a))
                trl60_misc5 = [trl60_misc5;bslsmp btsmp -hdr.Fs*2.5];
                
            % if misc4 tagged at 67
            elseif find(trig67_misc4 == allval(a))
                trl67_misc4 = [trl67_misc4;bslsmp btsmp -hdr.Fs*2.5];
                
            % if misc5 tagged at 67
            elseif find(trig67_misc5 == allval(a))
                trl67_misc5 = [trl67_misc5;bslsmp btsmp -hdr.Fs*2.5];
            end
        end
        a = a + 3;
        clear bslsmp sdsmp btsmp
    else 
        a = a + 1;
    end
    
end


