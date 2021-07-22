%% VS + RFT
% PhD project 2

% merge edf, fif and mat files 

% [c] Katharina Duecker
clear all; close all; clc; beep off
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
dtpth = fullfile(pth,'data');
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% document: subject id, n trial meg, el adjusted, response adjusted
mergesubj = cell(length(subjfolds),4);

for s = 1:length(subjfolds)
    
    mergesubj{s,1} = s;
    
    % load edf
    try
        d = dir(fullfile(dtpth,subjfolds{s}));
        f = {d.name};
        % find edf file
        idx = cellfun(@(x) regexp(x,'edf'),f,'UniformOutput',false);
        idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
        el = Edf2Mat(fullfile(dtpth,subjfolds{s},f{idxx}));
    catch ME
        mergesubj{s,3} = 'no el file';
        continue
    end
    
    % load events from maxfiltered data
    mgfls = dir(fullfile(pth,'results','meg', '1 maxfilter', subjfolds{s}));
    f = {mgfls.name};
    
    % find fif file
    idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    f = f(idxx);
    
    % load events separated in parts
    for fl = 1:length(f)
        events{fl} = ft_read_event(fullfile(pth,'results','meg', '1 maxfilter', subjfolds{s},f{fl}));
    end
    
    load(fullfile(pth,'experiment','trigdef.mat'))
    
    %% List trial type fif
    alltrl_bl = {};
    % trl: type, start sample, end sample
    for fl = 1:length(events)
        allval = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).value];   % these are all values, extract button press
        allsmp = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).sample];
        
        trl = [];
        a = 1;
        while a+2 <= length(allval)
            % a: beginning of sample, a+1: onset search display, a+2: button
            if allval(a) >= 5 && allval(a) <= trigdef{end,1} && allval(a+1) == 2
                
                if allval(a+2) > trigdef{end,1}
                    trl = [trl; allval(a) allsmp(a) allsmp(a+2)];
                else
                    trl = [trl; allval(a) allsmp(a) allsmp(a+1)+4*fs];
                end
            end
            a = a +1;
        end
        alltrl_bl{fl} = trl;
    end
    
    alltrl_list = vertcat(alltrl_bl{:});
    
    %% List trial type EL
    % extract el triggers
    idx = cellfun(@(x) regexp(x,'trig:'),el.Events.Messages.info,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    eltrials = el.Events.Messages.info(idxx);
    elsamples = el.Events.Messages.time(idxx);
    
    % extract trigger values
    extrtrig = cell2mat(cellfun(@str2num,cellfun(@(x) x(7:end),eltrials,'UniformOutput',false),'UniformOutput',false));
    
    % if practice trials are included (more than 24 1's = start of block) ->
    % only select experimental trials (last 24)
    
    begexp = find(extrtrig == 1);
    begexp = begexp(end-23:end);
    extrtrig = extrtrig(begexp(1):end);
    
    % eltrl: trigger, sample begin, sample end
    eltrl = [];
    a = 1;
    while a+2 <= length(extrtrig)
        % a: beginning of sample, a+1: onset search display, a+2: trial end (4)
        if extrtrig(a) >= 5 && extrtrig(a) <= trigdef{end,1} && extrtrig(a+1) == 2 && extrtrig(a+2) == 4
            eltrl = [eltrl;extrtrig(a) elsamples(a) elsamples(a+2)];
        end
        a = a + 1;
    end
    
    
    %% load mat file extract trial type and compare with edf & fif
    
    % 1. delete practice mat file
    
    f = {d.name};
    % find edf file
    idx = cellfun(@(x) regexp(x,'mat'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    f = f(idxx);
    if length(f) > 1
        for h = 1:length(f)
            if strfind(f{h},'pract')
                delete(fullfile(dtpth,subjfolds{s},f{h}))
            end
        end
    end
    
    
    load(fullfile(dtpth,subjfolds{s},f{1}))
    
    % 2. find specs & trigger!
    
    trltype = subj.exp.trials{1};
    
    rsptrig = zeros(size(trltype));
    
    % loop over block
    for b = 1:size(trltype,1)
        % loop over trials
        for t = 1:size(trltype,2)
            curspex = cellfun(@num2str,trltype{b,t},'UniformOutput',false);
            curspex = [curspex{:}];
            
            rsptrig(b,t) = trigdef{ismember(trigdef(:,2),curspex),1};
        end
    end
    
    % transpose for linear indexing
    rsptrig = rsptrig';
    
    %% Find common trials in all 3 - lead by MEG data
    
    % copy response and eyelink
    rsptrigc = reshape(rsptrig,[],1);
    eltrlc = eltrl(:,1);
    
    % which trial to reject?
    rejtrl_resp = [];                   % responses
    rejtrl_el = [];                     % eyelink
    rejtrl_meg = [];                    % meg
    
    w = 1;          % meg counter
    while w <= size(alltrl_list,1)
        if rsptrigc(w) == alltrl_list(w,1) && alltrl_list(w,1) == eltrlc(w)
            
        elseif rsptrigc(w) ~= alltrl_list(w,1) && alltrl_list(w,1) == eltrlc(w)
            % find first one that is the same
            fw = find(rsptrigc(w:end) == alltrl_list(w,1),1);
            rejtrl_resp = [rejtrl_resp, w:w+fw-2];
            rsptrigc(w:w+fw-2) = [];
            
        elseif rsptrigc(w) == alltrl_list(w,1) && alltrl_list(w,1) ~= eltrlc(w)
            % find first one that is the same
            fw = find(eltrlc(w:end) == alltrl_list(w,1),1);
            rejtrl_el = [rejtrl_el, w:w+fw-2];
            eltrlc(w:w+fw-2) = [];
            
        elseif rsptrigc(w) ~= alltrl_list(w,1) && alltrl_list(w,1) ~= eltrlc(w)
            fw = find(rsptrigc(w:end) == alltrl_list(w,1),1);
            rejtrl_resp = [rejtrl_resp, w:w+fw-2];
            rsptrigc(w:w+fw-2) = [];
            
            fw = find(eltrlc(w:end) == alltrl_list(w,1),1);
            rejtrl_el = [rejtrl_el, w:w+fw-2];
            eltrlc(w:w+fw-2) = [];
            
        else
            error('diff')
        end
        w = w + 1;
    end
    
    % create response cell
    rsp = cell(size(rsptrigc,1),3);
    % convert to cell
    rsp(:,1) = arrayfun(@(x) x*1,rsptrigc,'UniformOutput',false);
    
    % helper cell
    h = cell(prod(size(rsptrig)),2);
    % fill in RT and performance
    % loop over blocks
    c = 0;
    for b = 1:size(subj.rspns,2)
        for t = 1:size(subj.rspns{b},2)
            c = c + 1;
            h(c,:) = {subj.rspns{b}{end,t}, subj.rspns{b}{3,t}};
        end
    end
    
    h(rejtrl_resp,:) = [];
    
    rsp(:,2:end) = h;
    
    mergesubj{s,2} = length(alltrl_list);
    
    if length(eltrlc) ~= length(alltrl_list) || sum(abs(alltrl_list(:,1) - eltrlc)) > 0
        mergesubj{s,3} = 'check el';
    else
        mergesubj{s,3} = 'adjusted';
    end
    
    if length(rsptrigc) ~= length(alltrl_list) || sum(abs(alltrl_list(:,1) - rsptrigc)) > 0
        mergesubj{s,4} = 'check rsp';
    else
        mergesubj{s,4} = 'adjusted';
    end
    
    %% Find saccades
    
    % epoch EL data
    eltrl(rejtrl_el,:) = [];
    eltrl = [eltrl,zeros(size(eltrl,1),2)];
    elrt = zeros(size(eltrl,1),1);
    
    % xy coordinates over time: trials x time series position x x&y
    elxy = zeros(size(eltrl,1),5.5*1000,2);
    for e = 1:size(eltrl,1)
        % start sample
        bs = find(el.Samples.time == eltrl(e,2));
        es = find(el.Samples.time == eltrl(e,3));
        % store sample info for trial
        eltrl(e,end-1:end) = [bs es];
        elxy(e,1:es-bs+1,1) = el.Samples.posX(bs:es)';
        elxy(e,1:es-bs+1,2) = el.Samples.posY(bs:es)';

        % check rt
        elrt(e) = (es-bs)/1000 - 1.5;
        
        
    end
    
    % how many of the identified saccades and blinks would we have to
    % reject?
    
    % 1. saccades
    sacrej = [];
    for sac = 1:length(el.Events.Esacc.start)
        x = find((eltrl(:,4) <= el.Events.Esacc.start(sac)) + (eltrl(:,5) >= el.Events.Esacc.start(sac)) == 2 );
        y = find((eltrl(:,4) <= el.Events.Esacc.end(sac))  + (eltrl(:,5) >= el.Events.Esacc.end(sac))  == 2 );

        sacrej = [sacrej;x];
        sacrej = [sacrej;y];
    end
    sacrej = unique(sacrej);
    
    % blinks
    blinkrej = [];
    for sac = 1:length(el.Events.Eblink.start)
        x = find((eltrl(:,4) <= el.Events.Eblink.start(sac)) + (eltrl(:,5) >= el.Events.Eblink.start(sac)) == 2 );
        y = find((eltrl(:,4) <= el.Events.Eblink.end(sac))  + (eltrl(:,5) >= el.Events.Eblink.end(sac))  == 2 );
        
        blinkrej = [blinkrej;x];
        blinkrej = [blinkrej;y];
    end
    blinkrej = unique(blinkrej);
    
%     % check difference between el rt and response rt (~100 ms)
%     rtdiff = elrt-[rsp{:,3}]';
%     medrtdiff = median(rtdiff);
%     % relatively consistent difference
%     stdrtdiff = std(rtdiff);
    
%     % check x,y coord over time
%     plot(linspace(-1.5,4,1000*5.5),squeeze(elxy(1,:,1)))
%     hold on
%     plot(linspace(-1.5,4,1000*5.5),squeeze(elxy(1,:,2)))
%     legend('x','y')

    % store all el info in one struct
    elinfo.eltrl = eltrl;                   % trial info
    elinfo.el_rejtrl = rejtrl_el;           % reject based on meg
    elinfo.move = elxy;                     % movement over time
    elinfo.sac = sacrej;                    % identified trials with saccades
    elinfo.blink = blinkrej;                % identified trials with blinks
    
    % rsp info as struct
    rspinfo.trl = rsp;                      % trial info (trigger, performance, rt
    rspinfo.rsp_rejtrl = rejtrl_resp;       % rejected based on meg
    mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))
    save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
        'trl_overlap_meg_el_rsp.mat'),'trl','rspinfo','elinfo')
    
    
    clear events all* begexp a b d f* s t w el* rej* rsp* tr* idx* extrtrig curspex 
end
% meg file: [begin sample, end sample, baseline]


save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'),'mergesubj')
