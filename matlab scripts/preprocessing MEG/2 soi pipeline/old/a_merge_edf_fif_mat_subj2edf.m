%% VS + RFT
% PhD project 2

% merge edf, fif and mat files

% [c] Katharina Duecker
clear all; close all; clc; beep off
% define paths
pth = 'E:\UoB\Proj2_Visual Search';
mtpth = fullfile('Y:\Visual Search RFT','matlab scripts');
addpath(fullfile(mtpth,'edf-converter-master')); %edf2mat converter
% addpath('Y:\Visual Search RFT\matlab scripts\edfconverter')
addpath('C:\Program Files (x86)\SR Research\EyeLink\bin')
command = 'edf2asc.exe';
dtpth = fullfile(pth,'data');
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
mkdir(trl_merge_pth)
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% document: subject id, n trial meg, el adjusted, response adjusted
mergesubj = cell(length(subjfolds),4);

% screen settings: pixel -> degree
propixx_res = [1920 1080];                                   % propixx resolution
scr.w            = 72;                                       % screen width in cm
scr.h            = 40.5;                                     % screen height in cm
scr.d            = 142;
scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pi

%num_rej_trl_el = zeros(size(subjfolds));                           % number of trials contaminated by eye movements
for s = 20
    mergesubj{s,1} = subjfolds{s};
    sac_thresh = 1;
    
    try
        d = dir(fullfile(dtpth,subjfolds{s}));
        f = {d.name};
        % find edf file
        idx = cellfun(@(x) regexp(x,'.edf'),f,'UniformOutput',false);
        idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
        idxx = find(idxx);

        for fl = 1:length(f(idxx))
            filename = f{idxx(fl)};
            if fl == 1 & ~exist(fullfile(dtpth,subjfolds{s},[filename(1:end-4),'.asc']))
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
            end
            
            el_all{fl} = Edf2Mat(fullfile(dtpth,subjfolds{s},filename));
        end
    catch ME
        mergesubj{s,3} = 'no el file';
        continue
    end
    
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
    % List trial type
    % extract el triggers
    
    % create new structure el concatenating existing ones
    el = [];
    el.timeline = [el_all{1}.timeline;el_all{2}.timeline];
    el.normalizedTimeline = [el_all{1}.normalizedTimeline;el_all{2}.normalizedTimeline];
    el.Events.Messages.info = [el_all{1}.Events.Messages.info,el_all{2}.Events.Messages.info];
    el.Events.Messages.time = [el_all{1}.Events.Messages.time,el_all{1}.Events.Messages.time(end)+el_all{2}.Events.Messages.time];
    idx = cellfun(@(x) regexp(x,'trig:'),el.Events.Messages.info,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    if isempty(find(idxx))
        mergesubj{s,3} = 'no el file';
        clear events all* begexp a b d f fl t w el* rej* rsp* tr* idx* extrtrig curspex
        continue
    else
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
        a = 1;
        while a+2 <= length(extrtrig)
            % a: beginning of sample, a+1: onset search display, a+2: trial end (4)
            
            if extrtrig(a) >= 5 && extrtrig(a) <= trigdef{end,1} && extrtrig(a+1) == 2 && extrtrig(a+2) == 4
                eltrl = [eltrl;extrtrig(a) elsamples(a) elsamples(a+1) elsamples(a+2)];
            end
            a = a + 1;
        end
        
        
        % concatenate relevant el fields
        el.Samples.time = [el_all{1}.Samples.time;el_all{1}.Samples.time(end)+el_all{2}.Samples.time(end)];
        el.Samples.posX = [el_all{1}.Samples.posX; el_all{2}.Samples.posX];
        el.Samples.posY = [el_all{1}.Samples.posX; el_all{2}.Samples.posY];
        
        
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
                
                if isempty(subj.rspns{b}{3,t})
                    h(c,:) = {subj.rspns{b}{end,t}, 4};
                end
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
        
        % store all el info in one struct
        elinfo.eltrl = eltrl;                   % trial info
        elinfo.el_rejtrl = rejtrl_el;           % reject based on meg
        elinfo.move = elxy;                     % movement over time
       
        elinfo.num_rej_trl = size(unique([sacrej_liberal;blinkrej]),1);
        elinfo.num_rej_trl_lib = size(unique([sacrej_liberal;blinkrej]),1);

        elinfo.fsample = el_fsample;
        % rsp info as struct
        rspinfo.trl = rsp;                      % trial info (trigger, performance, rt
        rspinfo.rsp_rejtrl = rejtrl_resp;       % rejected based on meg
        mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))
        save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
            'trl_overlap_meg_el_rsp_sac_search.mat'),'alltrl_list','alltrl_bl','rspinfo','elinfo','-v7.3')
        
        
        clear events all* begexp a b d f fl t w el* rej* rsp* tr* idx* extrtrig curspex
    end
end
% meg file: [begin sample, end sample, baseline]

% num_rej_trl_el = zeros(size(subjfolds));  
% num_rej_trl_el_lib = zeros(size(subjfolds));  
% for s = 1:length(subjfolds)
%     if exist(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
%             'trl_overlap_meg_el_rsp_sac_search.mat'))
%         load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
%             'trl_overlap_meg_el_rsp_sac_search.mat'),'elinfo')
%         num_rej_trl_el(s) = elinfo.num_rej_trl;
%         num_rej_trl_el_lib(s) = elinfo.num_rej_trl_lib;
%     clear elinfo
%     end
%     
% end
% 
% save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_sac_search.mat'),'mergesubj','num_rej_trl_el','num_rej_trl_el_lib')
