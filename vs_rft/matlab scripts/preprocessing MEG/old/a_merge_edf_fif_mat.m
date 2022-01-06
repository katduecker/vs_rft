%% VS + RFT
% PhD project 2
% [c] Katharina Duecker

% merge edf, fif and mat files
% - script adjusts trial indices for the 3 files
% - if trial indices differ, user has to adjust trials by hand (too avoid
% errors and identify nature of inconsistencies)

% check results by follow-up script a2_check_merge

clear all; close all; clc; beep off

% define paths
pth = 'E:\UoB\Proj2_Visual Search';                                 % main path
mtpth = fullfile('Z:\Visual Search RFT','matlab scripts');          % matlab scripts path
addpath('C:\Users\katha\Documents\MATLAB\edf-converter-master');                    %edf2mat converter

addpath('C:\Program Files (x86)\SR Research\EyeLink\bin')           % eyelink path
command = 'edf2asc.exe';                                            % eyelink function edf2asc

dtpth = fullfile(pth,'data');                                       % data path
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');  % path where results are tb stored
mkdir(trl_merge_pth)

addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')                % fieldtrip
ft_defaults;

% list subjects
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% document: subject id, n trial meg, el adjusted, response adjusted
if exist(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))             % subject documentation exists already
    load(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))             % load
else
    mergesubj = cell(length(subjfolds),4);                                               % else create empty cell
end

noisy_subj = table2cell(readtable(fullfile(pth,'results','meg','1 maxfilter','noisy_subj.xlsx'),'ReadVariableNames',false));

propixx_res = [1920 1080];                                   % propixx resolution

btn = [5888, 6144];

for s = 3%:length(subjfolds)
    
    disp(['subj ', subjfolds{s}])
    el_fsample = 1000;
    mergesubj{s,1} = subjfolds{s};
    trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
    
    if ismember(noisy_subj,subjfolds{s})
        mergesubj{s,2} = 'MEG noisy';
    end

    % skip if done for this participant
%     if exist(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
%             'trl_overlap_meg_el_rsp.mat'))
%         disp('merge file exists')
%         continue
%     end

    %% Prep el file
    % if edf file has been converted to mat file already, load it
    if exist(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'))
        load(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'))
        
        
    else
        % convert edf file to mat struct using edf2mat
        try
            d = dir(fullfile(dtpth,subjfolds{s}));
            f = {d.name};
            % find edf file
            idx = cellfun(@(x) regexp(x,'.edf'),f,'UniformOutput',false);
            idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
            
            % if there is one edf file (experiment not interrupted)
            if length(find(idxx)) == 1
                if ~exist(fullfile(dtpth,subjfolds{s},[f{idxx}(1:end-4),'.asc']))
                    cd(fullfile(dtpth,subjfolds{s}))
                    system([command ' ' f{idxx}])
                    cd(mtpth)
                end
                el = Edf2Mat(fullfile(dtpth,subjfolds{s},f{idxx}));
                
                % use fieldtrip for sampling rate (this is very
                % inefficient...)
                
                % convert to asc
                cfg = [];
                cfg.dataset = fullfile(dtpth,subjfolds{s},[f{idxx}(1:end-4),'.asc']);
                el_ft = ft_preprocessing(cfg);
                
                
                if el_ft.fsample > 900
                    el_fsample = round(el_ft.fsample, -3);
                else
                    el_fsample = round(el_ft.fsample, -2);
                end
                clear el_ft
            else
                % more than one edf file -> experiment was interrupted at
                % some point
                
                idxx = find(idxx);
                for fl = 1:length(idxx)
                    filename = f{idxx(fl)};
                    if ~exist(fullfile(dtpth,subjfolds{s},[filename(1:end-4),'.asc']))
                        try
                            cd(fullfile(dtpth,subjfolds{s}))
                            system([command ' ' filename])
                            cd(mtpth)
                            
                            cfg = [];
                            cfg.dataset = fullfile(dtpth,subjfolds{s},[filename(1:end-4),'.asc']);
                            el_ft = ft_preprocessing(cfg);
                            if el_ft.fsample > 900
                                el_fsample = round(el_ft.fsample, -3);
                            else
                                el_fsample = round(el_ft.fsample, -2);
                            end
                            clear el_ft
                        catch ME
                        end
                    end
                    % convert all files to mat structures
                    el_all{fl} = Edf2Mat(fullfile(dtpth,subjfolds{s},filename));
                end
                
                % concatenate relevant el fields
                el = [];
                el.timeline = [el_all{1}.timeline;el_all{2}.timeline];
                el.normalizedTimeline = [el_all{1}.normalizedTimeline;el_all{2}.normalizedTimeline];
                el.Events.Messages.info = [el_all{1}.Events.Messages.info,el_all{2}.Events.Messages.info];
                el.Events.Messages.time = [el_all{1}.Events.Messages.time,el_all{1}.Events.Messages.time(end)+el_all{2}.Events.Messages.time];
                el.Samples.time = [el_all{1}.Samples.time;el_all{1}.Samples.time(end)+el_all{2}.Samples.time];
                el.Samples.posX = [el_all{1}.Samples.posX; el_all{2}.Samples.posX];
                el.Samples.posY = [el_all{1}.Samples.posX; el_all{2}.Samples.posY];
                el.Events.Esacc.start = [el_all{1}.Events.Esacc.start, el_all{2}.Events.Esacc.start];
                el.Events.Ssacc.time = [el_all{1}.Events.Ssacc.time,el_all{1}.Events.Ssacc.time(end)+el_all{2}.Events.Ssacc.time];
                el.Events.Esacc.duration = [el_all{1}.Events.Esacc.duration,el_all{2}.Events.Esacc.duration];
                el.Events.Sblink.time = [el_all{1}.Events.Sblink.time,el_all{2}.Events.Sblink.time];
            end
            
            % resample by hand if sampling rate higher than 1000 has been
            % used (by accident, this shouldn't happen)
            if el_fsample > 1000
                sr = el_fsample/1000;
                el_resamp = [];
                el_resamp.Samples.time = el.Samples.time(1:sr:length(el.Samples.time));
                el_resamp.Samples.posX = el.Samples.posX(1:sr:length(el.Samples.time));
                el_resamp.Samples.posY = el.Samples.posY(1:sr:length(el.Samples.time));
                el_resamp.Events.Esacc.start = el.Events.Esacc.start;
                el_resamp.Events.Ssacc.time = el.Events.Ssacc.time;
                el_resamp.Events.Esacc.duration  = el.Events.Esacc.duration;
                el_resamp.Events.Sblink.time = el.Events.Sblink.time; 
                el_resamp.Events.Messages.info = el.Events.Messages.info;
                el_resamp.Events.Messages.time = el.Events.Messages.time;
                
                el = el_resamp;
                save(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'),'el','-v7.3')

            end
            if ~exist(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'))
                save(fullfile(dtpth,subjfolds{s},'el_edf2mat.mat'),'el','-v7.3')
            end
        catch ME
            disp('no el file');
            continue
        end
    end
    
    % if elinfo already exists (re-running script, skip above disabled)
    % load el to skip all of the above
    if exist(fullfile(trl_merge_pth,subjfolds{s},'trl_overlap_meg_el_rsp_sac_search.mat'))
        load(fullfile(trl_merge_pth,subjfolds{s},'trl_overlap_meg_el_rsp_sac_search.mat'),'elinfo')
        
        try
            el_fsample = elinfo.fsample;
            clear elinfo
        catch ME
        end
    end
    
    %% Prep MEG 
    % load events from maxfiltered data
    mgfls = dir(fullfile(dtpth,subjfolds{s},'meg'));
    f = {mgfls.name};
    
    % find fif file
    idx = cellfun(@(x) regexp(x,'part'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    f = f(idxx);
    
    % load events separated in parts
    for fl = 1:length(f)
        events{fl} = ft_read_event(fullfile(dtpth,subjfolds{s},'meg',f{fl}));
    end
    
    load(fullfile(pth,'experiment','trigdef.mat'))
    
    %  List trial type fif
    alltrl_bl = {};
    % trl: type, start sample, end sample
    for fl = 1:length(events)
        allval = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).value];   % these are all values, extract button press
        allsmp = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).sample];
        
        % get trigger indices
        %trig_bsl = find((allval >= 5)+(allval<=trigdef{end,1}) == 2);
        trig_bsln= ((allval >= 5)+(allval<=trigdef{end,1}) == 2);
        trig_bslb1 = ((allval-btn(1) >= 5)+(allval-btn(1)<=trigdef{end,1}) == 2);
        trig_bslb2 = ((allval-btn(2) >= 5)+(allval-btn(2)<=trigdef{end,1}) == 2);

        trig_bsl = find(logical(trig_bsln + trig_bslb1 + trig_bslb2));
        find((trig_bsln + trig_bslb1 + trig_bslb2)==2)
%         % check triggers that are not followed by 2 (begin of search
%         % display)
         % if there is 2 trial start trigger behind each other (why?!)
         % delete
          notrl_idx = (allval(trig_bsl+1) >= 5)+(allval(trig_bsl+1)<=trigdef{end,1}) == 2;
%         norsp_idx = find((allval(trig_bsl+2) <2000) + (allval(trig_bsl+2)~=4) == 2);
         trig_bsl(unique(notrl_idx)) = [];

        trl = [allval(trig_bsl)',allsmp(trig_bsl)' allsmp(trig_bsl+1)' allsmp(trig_bsl+2)'];
        trigtrl{fl} = [allval(trig_bsl)',allval(trig_bsl+1)',allval(trig_bsl+2)',allval(trig_bsl+3)']
        alltrl_bl{fl} = trl;
    end
    
    % if there is not enough identified trials, chances are the triggers
    % were messed up by P pressing button for too long
    alltrl_list = vertcat(alltrl_bl{:});
    trigtrl_list = vertcat(trigtrl{:});
    if length(alltrl_list) ~= 960
        mergesubj{s,3} = 'trigger messed up';
        continue
    end

    meg_rt = (alltrl_list(:,4) - alltrl_list(:,3))/1000;
    
    %% List trial type EL
    % extract el triggers

    idx = cellfun(@(x) regexp(x,'trig:'),el.Events.Messages.info,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    if isempty(find(idxx))
        mergesubj{s,3} = 'no el file';
        clear events all* begexp a b d f fl t w el* rej* rsp* tr* idx* extrtrig curspex
        continue
    end
    eltrials = el.Events.Messages.info(idxx);
    elsamples = el.Events.Messages.time(idxx);
    
    % extract trigger values
    extrtrig = cell2mat(cellfun(@str2num,cellfun(@(x) x(7:end),eltrials,'UniformOutput',false),'UniformOutput',false));
    
    % if practice trials are included (more than 24 1's = start of block) ->
    % only select experimental trials (last 24)
    
    begexp = find(extrtrig == 1);
    begexp = begexp(end-23:end);
    extrtrig = extrtrig(begexp(1):end);
    elsamples = elsamples(begexp(1):end);
    eltrials = eltrials(begexp(1):end);
    
    % eltrl: trigger, sample begin, sample end
    
    
    eltrl = [];
    trig_bsl = 1;
    while trig_bsl+2 <= length(extrtrig)
        % a: beginning of sample, a+1: onset search display, a+2: trial end (4)
        
        if extrtrig(trig_bsl) >= 5 && extrtrig(trig_bsl) <= trigdef{end,1} && extrtrig(trig_bsl+1) == 2 && extrtrig(trig_bsl+2) == 4
            eltrl = [eltrl;extrtrig(trig_bsl) elsamples(trig_bsl) elsamples(trig_bsl+1) elsamples(trig_bsl+2)];
        end
        trig_bsl = trig_bsl + 1;
    end
    
    % reaction time according to eye tracker
    eltrl = [eltrl,zeros(size(eltrl,1),3)];
    elrt = zeros(size(eltrl,1),1);
    
    % extract samples corresponding to time
    fixdot = propixx_res./2; % coordinates fixation do
    elxy = zeros(size(eltrl,1),5.5.*el_fsample,2);
    for e = 1:size(eltrl,1)
        % start sample: index of time
        bs = find(el.Samples.time == eltrl(e,2));
        vs = find(el.Samples.time == eltrl(e,3));
        es = find(el.Samples.time == eltrl(e,4));
        % store sample info for trial
        eltrl(e,end-2:end) = [bs vs es];
        elxy(e,1:es-bs+1,1) = el.Samples.posX(bs:es)';
        elxy(e,1:es-bs+1,2) = el.Samples.posY(bs:es)';
        
        % check rt
        elrt(e) = (eltrl(e,4)-eltrl(e,3))/el_fsample;
 
    end
    
    %% Trial type mat file (responses stored during experiment)
    
    % delete practice mat file
    d = dir(fullfile(dtpth,subjfolds{s}));
    
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
    
    d = dir(fullfile(dtpth,subjfolds{s}));

    f = {d.name};
    idx = cellfun(@(x) regexp(x,'mat'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    f = f(idxx);
    
    % find response files
    idx = cellfun(@(x) regexp(x,subjfolds{s}(end-3:end)),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    f = f(idxx);
    
    
    % if there are 2 response files, concatenate
    if length(f) == 1
        load(fullfile(dtpth,subjfolds{s},f{1}))
    else
        for i = 1:length(f)
            subj_all{i} = load(fullfile(dtpth,subjfolds{s},f{i}));
        end
        subj = subj_all{1}.subj;
        
        subj.exp.trials{1} = [];
        subj.exp.trials{2} = [];
        subj.srchdsp = {};
        subj.rspns = {};
        for i = 1:length(f)
            subj.exp.trials{1} = [subj.exp.trials{1};subj_all{i}.subj.exp.trials{1}];
            subj.exp.trials{2} = [subj.exp.trials{2};subj_all{i}.subj.exp.trials{2}];
            subj.srchdsp = [subj.srchdsp;subj_all{i}.subj.srchdsp];
            subj.rspns = [subj.rspns,subj_all{i}.subj.rspns];
        end
        
        % number of blocks = 23
        subj.exp.trials{1} = subj.exp.trials{1}(end-23:end,:);
        subj.exp.trials{2} = subj.exp.trials{2}(end-23:end,:);
        subj.rspns = subj.rspns(end-23:end);
        subj.srchdsp = subj.srchdsp(end-23:end,:);
    end
    
    % find trial specifics & trigger!
    
    trltype = subj.exp.trials{1};
    
    rsptrig = zeros(size(trltype));
    rt = zeros(size(trltype));
    % loop over block
    for b = 1:size(trltype,1)
        % loop over trials
        for t = 1:size(trltype,2)
            curspex = cellfun(@num2str,trltype{b,t},'UniformOutput',false);
            curspex = [curspex{:}];
            
            rsptrig(b,t) = trigdef{ismember(trigdef(:,2),curspex),1};
            % store reaction time
            if ~isempty(subj.rspns{b}{3,t})
                rt(b,t) = subj.rspns{b}{3,t};
            else 
                rt(b,t) = 4;
            end
        end
    end
    
    % transpose for linear indexing
    rsptrig = rsptrig';
    rt = rt';
    
    rsptrigc = reshape(rsptrig,[],1);
    rtc = reshape(rt,[],1);
    
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
            
            if isempty(subj.rspns{b}{3,t})
                h(c,:) = {subj.rspns{b}{end,t}, 4};
            end
        end
    end
       
    rsp(:,2:end) = h;
    
    isequal([rsp{:,3}]',rtc) % sanity check
    isequal([rsp{:,1}]',rsptrigc) % sanity check

    % comparison rtc and meg_rt: only 1st two decimal places are the same
    rtc = floor(rtc.*100)./100;
    meg_rt = floor(meg_rt.*100)./100;
    
    %% Find common trials in all 3 - adjust by hand
    
    % copy response and eyelink

    eltrlc = eltrl(:,1);
    alltrl_listc = alltrl_list;
    
    % which trial to reject?
    rejtrl_resp = [];                   % responses
    rejtrl_el = [];                     % eyelink
    rejtrl_meg = [];                    % meg
    wm = 1;          % meg counter
    wr = 1;          % response counter      
    % This loop assumes that both the responses and the eyetracker contain
    % all trials! (which should be the case)
    while wm <= size(alltrl_listc,1)
        if rsptrigc(wr) ~= alltrl_list(wm,1) && alltrl_list(wm,1) ~= eltrlc(wr)
            % ideally show the next 10, if not possible, show until end

            % find index
            try
                x = find(alltrl_listc(wm,1) == rsptrigc(wr:wr+10),1);
            catch
                x = find(alltrl_listc(wm,1) == rsptrigc(wr:end),1);

            end
            % check if rt is also the same
            x_rt = meg_rt(wm) >= rtc(wr+x-1)-0.14 + meg_rt(wm) <= rtc(wr+x-1)+0.14;
            if x_rt
                fw_rt = x-1;
                fw_meg = 0;

            else
                if wm+10 > size(alltrl_listc,1)
                    [rsptrigc(wr:end),rtc(wr:end)]
                    [alltrl_listc(wm:end,1),meg_rt(wm:end)]

                else
                    [rsptrigc(wr:wr+10),alltrl_listc(wm:wm+10,1), rtc(wr:wr+10), meg_rt(wm:wm+10)]
                end
                fw_rt = input('How many to delete rsp?')
                fw_meg = input('How many to delete meg?')


            end


            rejtrl_resp = [rejtrl_resp, wr:wr+fw_rt-1];
            
            rejtrl_el = [rejtrl_el, wr:wr+fw_rt-1];
   
            % increase counter response
            wr = wr+fw_rt+1;
            wm = wm+fw_meg+1;              % increase counter meg
        else
            % increase counter
            wr = wr+1;
            wm = wm+1;

        end
    end
    %break
    
    
    %% save info in elinfo, meginfo, rspinfo
    
    mergesubj{s,2} = length(alltrl_listc);
    
    % store all el info in one struct
    elinfo.eltrl = eltrl;                   % trial info
    elinfo.el_rejtrl = rejtrl_el;           % reject based on meg
    elinfo.move = elxy;                     % movement over time
    meginfo.rejtrl_meg = rejtrl_meg;
    meginfo.meg_rt = meg_rt;
    meginfo.trigtrl = trigtrl_list;
    meginfo.alltrl_list  = alltrl_list;
    meginfo.alltrl_bl = alltrl_bl;
    elinfo.fsample = el_fsample;
    % rsp info as struct
    rspinfo.trl = rsp;                      % trial info (trigger, performance, rt
    rspinfo.rsp_rejtrl = rejtrl_resp;       % rejected based on meg
    
    % discard RT <200 ms & = 4 sec
    rtc(rejtrl_resp) = [];
    rspinfo.rejfast = find(rtc<=0.2);
    rspinfo.rejslow = find(rtc==4);
    mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))
    

    
    % check if all went okay!
    [trials_ok, rt_ok, idx_rt_diff] = a1_check_merged_trls(rspinfo, meginfo);
    % if it didn't, use a2_check_merge.m!

    if ~trials_ok || ~rt_ok
        mergesubj{s,3} = 'trials and RT mismatch';
    else
        mergesubj{s,3} = 'adjusted';
        mergesubj{s,4} = 'adjusted';
        
    end
    mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))
    save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
        'trl_overlap_meg_el_rsp.mat'),'rspinfo','elinfo','meginfo','idx_rt_diff','-v7.3')
    
   
    
   % clear events all* begexp a b d f fl t w el* rej* rsp* tr* idx* extrtrig curspex
    save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'),'mergesubj')

    
end

save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_sac_search.mat'),'mergesubj')
