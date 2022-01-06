%% VS + RFT
% PhD project 2
% [c] Katharina Duecker

% merge edf, fif and mat files
% - script adjusts trial indices for the 3 files
% - if trial indices differ, user has to adjust trials by hand (too avoid
% errors and identify nature of inconsistencies)

clear all; close all; clc; beep off

% define paths
pth = 'D:\UoB\Proj2_Visual Search';                                 % main path
mtpth = fullfile('Z:\Visual Search RFT','matlab scripts');          % matlab scripts path
% addpath('C:\Users\katha\Documents\MATLAB\edf-converter-master');                    %edf2mat converter
% 
% addpath('C:\Program Files (x86)\SR Research\EyeLink\bin')           % eyelink path
command = 'edf2asc.exe';                                            % eyelink function edf2asc

dtpth = fullfile(pth,'data');                                       % data path
maxfpth = fullfile(pth, 'results','meg','1 maxfilter');
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');  % path where results are tb stored
mkdir(trl_merge_pth)

addpath('C:\Users\dueckerk\Documents\MATLAB\fieldtrip\')                % fieldtrip
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

load(fullfile(pth,'experiment','button_val.mat'),'buttons_right');
buttons_right = buttons_right(1:3);
for s = 41%length(subjfolds)
    
    disp(['subj ', subjfolds{s}])
    el_fsample = 1000;
    mergesubj{s,1} = subjfolds{s};
    trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
    
    if ~isempty(find(ismember(noisy_subj,subjfolds{s})))
        mergesubj{s,2} = 'MEG noisy';
        continue
    end

    %skip if done for this participant
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
                save(fullfile(dtpth,subjfolds{s},'el_struct.mat'),'el','-v7.3')

            end

        catch ME
            mergesubj{s,3} = 'no el file';

            disp('no el file');
            continue
        end
    end
    
    % if elinfo already exists (re-running script, skip above disabled)
    % load el to skip all of the above
    if exist(fullfile(trl_merge_pth,subjfolds{s},'trl_overlap_meg_el_rsp.mat'))
        load(fullfile(trl_merge_pth,subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'elinfo')
        
        try
            el_fsample = elinfo.fsample;
            clear elinfo
        catch ME
        end
    end
    
    %% Prep MEG 
    % load events from maxfiltered data
    mgfls = dir(fullfile(maxfpth,subjfolds{s}));
    f = {mgfls.name};
    
    % find fif file
    idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    f = f(idxx);
    
    % load events separated in parts
    for fl = 1:length(f)
        events{fl} = ft_read_event(fullfile(maxfpth,subjfolds{s},f{fl}));
    end
    
    load(fullfile(pth,'experiment','trigdef.mat'))
    
    %  List trial type fif
    alltrl_bl = {};
    trigtrl = {};
    % trl: type, start sample, end sample
    for fl = 1:length(events)
        allval = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).value];   % these are all values, extract button press
        allsmp = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).sample];

        % get trigger indices
        trig_bsl = find((allval >= 5)+(allval<=trigdef{end,1}) == 2);

        % if there is 2 trial start trigger behind each other (why?!)
        % delete
        notrl_idx = (allval(trig_bsl+1) >= 5)+(allval(trig_bsl+1)<=trigdef{end,1}) == 2;
        trig_bsl = trig_bsl(~notrl_idx);

        trl = [allval(trig_bsl)',allsmp(trig_bsl)' allsmp(trig_bsl+1)' allsmp(trig_bsl+2)'];
        trigtrl{fl} = [allval(trig_bsl)',allval(trig_bsl+1)',allval(trig_bsl+2)'];
        alltrl_bl{fl} = trl;
    end
    
    % if there is not enough identified trials, chances are the triggers
    % were messed up by P pressing button for too long
    alltrl_list = vertcat(alltrl_bl{:});
    trigtrl_list = vertcat(trigtrl{:});
  

    if length(trigtrl_list) < 800
                mergesubj{s,3} = 'trigger messed up';
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

    rsptrig = rsptrig';   
    rsptrigc = reshape(rsptrig,[],1);
    
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
    
    meg_rt = floor(meg_rt.*100)./100;
    
    if length(rsp) > length(trigtrl_list)
        wm = 1;
        wr = 1;
        rej_rsp = [];
        rsp_trig = [rsp{:,1}]';
        while wm <= length(trigtrl_list) && wr < length(rsptrigc)
            if rsp_trig(wr) ~= trigtrl_list(wm,1)
                % find the first one that is the same
                x = find(rsp_trig(wr:end)==trigtrl_list(wm,1),1);
                rej_rsp = [rej_rsp,wr:wr+x-2];
                wr = wr + x -1;
            else
                wr = wr +1;
                wm = wm +1;
            end

        end
    elseif length(rsp) < length(trigtrl_list)
        error('MEG triggers weird')
    else
        rej_rsp = [];
    end

    %% Trial rejection

    % MEG
    % Reject MEG trials that contain button presses at wrong time (P kept pressing buttons)
    meg_rej = logical((trigtrl_list(:,1) > 36) + (trigtrl_list(:,2) ~=2));

    rt_array = [rsp{:,3}]';
    rt_array(rej_rsp) = [];
    % RT
    % reject that contain RT < 200 ms or RT > 4 s
    rt_rej = logical(rt_array< 0.2 + rt_array == 4);

    % mutual reject
    rej_trl = logical(meg_rej + rt_rej);

    % difference stim RT and MEG RT
    rt_diff = rt_array-meg_rt;
    idx_bigdiff = abs(rt_diff) > 0.014;
    rej_trl = logical(rej_trl+idx_bigdiff);

    %% save info in elinfo, meginfo, rspinfo
    
    mergesubj{s,2} = length(alltrl_list);
    
    % store all el info in one struct
    elinfo.eltrl = eltrl;                   % trial info
    elinfo.move = elxy;                     % movement over time
    elinfo.rejtrl_all = rej_trl;            % reject these for all
    elinfo.rej_rsp = rej_rsp;               % reject specifically for eyelink and rsp
    meginfo.rejtrl_all = rej_trl;           % reject for all
    meginfo.meg_rt = meg_rt;
    meginfo.trigtrl = trigtrl_list;
    meginfo.alltrl_list  = alltrl_list;
    meginfo.alltrl_bl = alltrl_bl;
    elinfo.fsample = el_fsample;
    % rsp info as struct
    rspinfo.trl = rsp;                      % trial info (trigger, performance, rt
    rspinfo.rejtrl_all = rej_trl;       % rejected for all
    rspinfo.rejtrl_rsp = rej_rsp;       % reject only for response

    mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))


    mergesubj{s,3} = 'adjusted';
    mergesubj{s,4} = 'adjusted';

    mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))
    save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
        'trl_overlap_meg_el_rsp.mat'),'rspinfo','elinfo','meginfo','-v7.3')
    
   
    
    save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'),'mergesubj')
    clear events all* begexp a b d f fl t w el* rej* rsp* tr* idx* extrtrig curspex

    
end

save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_sac_search.mat'),'mergesubj')
