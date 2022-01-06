%%%% 2020-02-23
%%%% using hilbert to compute coherence
%%%% for combined dataset


function TagCoh_hilbert(sid)
%%% choosing significant sensors based on epoch type eee
%%% using epoch aligned with rrr
%%% running dataset ddd
tic
server = 1;
eee = 1; %%% epoch index for sensor selection
rrr = 1; %%% running condition index
ddd = 1; %%% dataset condition index; for cmb-datasete, just use the name of Targ60

%%%% running setting
runconds = {'WrdOn','WrdOff'};
RunCond = runconds{rrr};
DataSets = {'Targ60','JEP60'};
EpochType = {'PreTarg';'Targ';'PosTarg'};
targid = [-1 0 1];
tagfreq = 60;
freqrange = 40:2:80;
freq_halfwidth = 5; %%% for hilbert filter
tagfreq_id = 11;
conds = {'all';'low';'high'};

%%%% basic settingup
if server
    rootdir = '/rds/projects/2018/jenseno-reading/';
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
    addpath(genpath('/rds/projects/2018/jenseno-reading/Analyse_codes/helper_functions/plotfig/'))
else
    %%% local
    rootdir = 'Z:\';
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    addpath(genpath('Z:\Analyse_codes\helper_functions\plotfig\'));
end
%%% result path
PPath.FigPath = [rootdir 'Results' filesep 'Coh_hilbert_01' filesep 'minirt_trl' filesep];
% PPath.FigPath = [rootdir 'Results' filesep 'WrdOff' filesep];

%%% get the subinformation
load([rootdir 'Analyse_data' filesep 'SubInfo.mat']);
sss_targ = SubInfo.sub_has2sets_id_inTarg60;
%%%+++ get the subs that have 2 datasets
% % % % % % subtarg = SubInfo.subjects.Targ60;
% % % % % % subjep = SubInfo.subjects.JEP60;
% % % % % % sss_targ = []; sss_jep = [];
% % % % % % for tmpsub = 1:length(subjep)
% % % % % %     stmp = find(strcmp(subjep{tmpsub}(1:end-4),subtarg));
% % % % % %     if stmp
% % % % % %         sss_targ = [sss_targ; stmp];
% % % % % %         sss_jep = [sss_jep; tmpsub];
% % % % % %     end
% % % % % % end
% % % % % % %%% save sub_idx that have both datasets
% % % % % % SubInfo.sub_has2sets_id_inTarg60 = sss_targ;
% % % % % % save([rootdir 'Analyse_data' filesep 'SubInfo.mat'],'SubInfo');


if server == 1
    %%% set output
    TagCoh = [];
    TagCoh.conds = conds;
    TagCoh.freq_halfwidth = freq_halfwidth;
    
    %%%%%%% defining datasets
    eval(['subjects = SubInfo.subjects.' DataSets{ddd} ';']);
    
    %%%%%%===================== dataset loop
    TagCoh.subs = subjects;
    TagCoh.SigTagSen_Occip = cell(size(subjects));
    TagCoh.Coh_hdr = {'freq*time*sub*cond'};
    
    %%%%%%===================== subject loop
    for sss = sss_targ(sid) %% 1:length(subjects)
        fprintf(['*******************analyzing: ' DataSets{ddd} '-s' num2str(sss) '*********************** \n\n']);
        sub = subjects{sss};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%% 1. set all kinds of paths and filenames, loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%======== JEP 60
        PPath.SaveData = [rootdir 'Analyse_data' filesep sub '_JEP' filesep];
        PPath.SubPath = [PPath.SaveData DataSets{2} '_TW1000' filesep];
        load([PPath.SubPath 'epoch_' RunCond]);
        eval(['epochdata_jep = epoch_' RunCond ';']);
        %%%==== targ 60
        PPath.SaveData = [rootdir 'Analyse_data' filesep sub filesep];
        PPath.SubPath = [PPath.SaveData DataSets{1} '_TW1000' filesep];
        load([PPath.SubPath 'epoch_' RunCond]);
        eval(['epochdata = epoch_' RunCond ';']);
        %%%=== combining 2 datasets together
        epochdata.trialinfo = [epochdata.trialinfo; epochdata_jep.trialinfo];
        epochdata.trial = [epochdata.trial epochdata_jep.trial];
        epochdata.time = [epochdata.time epochdata_jep.time];
        clear epoch_* epochdata_jep
        
        %%%%===== normalize pd and remove trials that without photodies signal,
        %%%% no saving out, rawdata is rawdata
        %%% get pd label index
        pdi = [find(strcmp(epochdata.label,'MISC004')) find(strcmp(epochdata.label,'MISC005'))];
        rmtrl = [];
        for ttt = 1:length(epochdata.trial)
            %%% remove trials that have no pd-004 signal
            if max(epochdata.trial{ttt}((pdi(1)),:)) < 0.005
                rmtrl = [rmtrl; ttt];
            end
            epochdata.trial{ttt}((pdi),:) = zscore(epochdata.trial{ttt}((pdi),:));
        end
        %%% get first fixation trials
        %%% event_raw_header = {'sentence_id','word_loc','loc2targ','word_freq','word_length','saccade2this_duration','fixation_on_MEG','fixation_duration','NextOrder','FirstPassFix','PreviousOrder','SentenceCondition'};
        %%% NextOrder:location of next fixated word compared to current word; FirstPassFix: fixate in the 1st pass time
        validduration = find(epochdata.trialinfo(:,10) == 1); %% epochdata.trialinfo(:,9) >=0 &
        cfg = [];
        cfg.trials = setdiff(validduration,rmtrl);
        epochdata = ft_selectdata(cfg, epochdata);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. getting significant tagging response sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% permutation test between epoch data and baseline data
        %%% selecting comparing epochs
        flickid = targid(eee);
        cfg = [];
        cfg.trials = find(epochdata.trialinfo(:,3) == flickid);
        data_all = ft_selectdata(cfg, epochdata);
        load([PPath.SubPath 'epoch_BL_Cross']);
        
        %%%======= get fourier spctrm --- no time domain
        cfg              = [];
        cfg.output       = 'fourier';
        cfg.channel      = {'MEGGRAD','MISC004'};
        cfg.method       = 'mtmfft';
        cfg.taper        = 'hanning';
        cfg.foi          = tagfreq;
        cfg.pad          = 'nextpow2';
        fourier          = ft_freqanalysis(cfg,data_all);%=% tfr.powspctrm:3D data: trialnum*channelnum*freqpoint
        fourier_bl       = ft_freqanalysis(cfg,epoch_BL_Cross);
        
        %%%% statistical test of coherence
        %%%% https://mailman.science.ru.nl/pipermail/fieldtrip/2018-January/012000.html
        %%%% --- doesnot work! can not seperate the
        %%%% channelcmb in stat,it's huge:
        %%%% since 'cluster correction' is not suitable here in the coherence
        %%%% statistic,(https://mailman.science.ru.nl/pipermail/fieldtrip/2010-April/002830.html)
        %%%% and the stat cann't output by individual channel combination, I trick it as single channel combination each time manually
        load([rootdir 'OccipSens_full'])
        nchancom = length(OccipSens_full);
        stat_mask = zeros(nchancom,1);
        
        tic
        for i = 1:nchancom
            cfg            = [];
            cfg.channel    = {OccipSens_full{i},'MISC004'};
            fourier_tmp    = ft_selectdata(cfg, fourier);
            fourier_bl_tmp = ft_selectdata(cfg, fourier_bl);
            
            cfg                  = [];
            cfg.parameter        = 'fourierspctrm';
            cfg.frequency        = tagfreq;
            cfg.statistic        = 'ft_statfun_indepsamplesZcoh';  %%%% take fourierspctrm as input, so no time domain information
            cfg.method           = 'montecarlo';
            cfg.tail             = 1; %% right sided, grp1 is bigger than grp2
            cfg.alpha            = 0.01;
            cfg.numrandomization = 10000;
            ntrl_1 = size(fourier.fourierspctrm,1);
            ntrl_2 = size(fourier_bl.fourierspctrm,1);
            design = zeros(1, ntrl_1 + ntrl_2);
            design(1,1:ntrl_1) = 1;
            design(1,(ntrl_1 +1):(ntrl_1 + ntrl_2))= 2;
            cfg.design = design;
            cfg.ivar   = 1;
            cfg.design = design;
            stat = ft_freqstatistics(cfg, fourier_tmp, fourier_bl_tmp);
            stat_mask(i) = stat.mask;
            stat_p(i,1)  = stat.prob;
        end
        toc
        TagCoh.Stat(:,1) = stat_mask;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. runing the 2D coherence in different conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TagCoh.freq = freqrange;
        TagCoh.time = epochdata.time{1};
        %%%===== epoch and conditon loop
        for mmm = 1:length(EpochType)
            %%% get equal trial number across all conditions
            %%%% seperate trials based on target-freq
            idx_low = find(epochdata.trialinfo(:,3) == targid(mmm) & epochdata.trialinfo(:,12) == 1);
            idx_high = find(epochdata.trialinfo(:,3) == targid(mmm) & epochdata.trialinfo(:,12) == 2);
            nlow = length(idx_low);
            nhigh = length(idx_high);
            if nlow > nhigh
                idx_low = idx_low(randperm(nlow));
                idx_low = idx_low(1:nhigh);
            elseif nlow < nhigh
                idx_high= idx_high(randperm(nhigh));
                idx_high = idx_high(1:nlow);
            end
            
            %%% get epochs of all conditions
            cfg        = [];
            cfg.trials = idx_low;
            data_low   = ft_selectdata(cfg, epochdata);
            cfg.trials = idx_high;
            data_high  = ft_selectdata(cfg, epochdata);
            cfg        = [];
            cfg.trials = find(epochdata.trialinfo(:,3) == targid(mmm)); %& epochdata.trialinfo(:,8) >= 4*1000/tagfreq);
            data_all   = ft_selectdata(cfg, epochdata);
            idx_all    = 1:size(data_all.trialinfo,1);
            
            %%% get trl numbers in each condition
            TrlNum = [size(data_all.trial,2) nlow nhigh];
            eval(['TagCoh.' EpochType{mmm} '_TrlNum_raw(1,:) = TrlNum;']);
            for ccc = 1:length(conds)
                eval(['TagCoh.' EpochType{mmm} '_TrlInfo{1,ccc} = data_' conds{ccc} '.trialinfo;']);
                eval(['TagCoh.' EpochType{mmm} '_RT(1,ccc) = nanmean(data_' conds{ccc} '.trialinfo(:,8));']);
                eval(['TagCoh.' EpochType{mmm} '_TrlIdx_equ{1,ccc} = idx_' conds{ccc} ';']);
            end
            
            %%%%%%%%%%%=================== Coherence===============%%%%%%%%%%%%
            %%%% conditon loop
            for ccc = 1:length(conds)
                %%%% just to get a dummy coherence struct
                cfg            = [];
                cfg.output     = 'fourier';
                cfg.channel    = {'MEGGRAD','MISC004'};
                cfg.method     = 'mtmconvol';
                cfg.taper      = 'hanning';
                cfg.foi        = freqrange;
                cfg.toi        = -0.5:0.5:0.5;
                cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
                cfg.keeptrials = 'yes';
                cfg.pad        = 'nextpow2';
                eval(['fourier = ft_freqanalysis(cfg, data_' conds{ccc} ');']);  
                %%% get coherence spctrm
                cfg            = [];
                cfg.method     = 'coh';
                cfg.channelcmb = {'MEGGRAD', 'MISC004'}; %% all channell combinations calculated together
                coh            = ft_connectivityanalysis(cfg, fourier); %% the raw coherence
                coh.label      = fourier.label(1:length(coh.labelcmb));
                coh.time       = epochdata.time{1};
                
                %%% do the real coherence using hilbert complex
                coh.cohspctrm = zeros(length(coh.label),length(coh.freq),length(coh.time));
                for fff = 1:length(freqrange)
                    %%%% get hilbert filert data
                    cfg = [];
                    cfg.channel    = {'MEGGRAD','MISC004'};
                    cfg.bpfilter   = 'yes';
                    cfg.bpfreq     = [freqrange(fff)-freq_halfwidth freqrange(fff)+freq_halfwidth];
                    cfg.hilbert    = 'complex';
                    cfg.keeptrials = 'yes';
                    eval(['fltdata = ft_preprocessing(cfg,data_' conds{ccc} ');']);
                    for chan = 1:length(fltdata.label)-1
                        for ttt = 1:length(fltdata.trial)
                            sig1(:,ttt) = fltdata.trial{ttt}(chan,:);
                            sig2(:,ttt) = fltdata.trial{ttt}(end,:);
                        end
                        spec1 = nanmean(sig1.*conj(sig1),2);
                        spec2 = nanmean(sig2.*conj(sig2),2);
                        specX = abs(nanmean(sig1.*conj(sig2),2)).^2;
                        coh.cohspctrm(chan,fff,:) = specX./(spec1.*spec2);
                    end
                end
                
                %%% finally get the sig occipital tagging sensors
                if mmm == 1 && ccc == 1
                    sigsens = OccipSens_full(stat_mask > 0);
                    %%% get all the sig tagging sensors in the order of p value
                    sigp   = stat_p(stat_mask > 0);
                    [~, Isig] = sort(sigp,'ascend');
                    TagCoh.SigTagSen_All{1,1} = sigsens(Isig);
                    %%%%% get the strongest responder
                    tmp = nanmean(nanmean(coh.cohspctrm(:,tagfreq_id,:),2),3);
                    %%% set nan to negtive value, otherwise nan would be the
                    %%% biggest value
                    tmp(isnan(tmp)) = 0;
                    TagCoh.SigTagSen_Occip(1,1) = {'No'};
                    TagCoh.SigTagSen_Occip{1,2} = 0;
                    sigchanid = 156;
                    [~, I] = sort(tmp,'descend');
                    sigsen = [];
                    %%%% first n strongest tagging sensors
                    for kkk = 1:length(I)  %%% check the strongest n sensors
                        tmpsen = coh.label(I(kkk)); %%% the strongest tagging responder
                        if ~isempty(sigsens) && ismember(tmpsen,sigsens)
                            sigsen = [sigsen, tmpsen];
                            sigchanid(1,length(sigsen)) = I(kkk);
                            TagCoh.SigTagSen_Occip{1,1} = sigsen;
                            TagCoh.SigTagSen_Occip{1,2} = sigchanid;
                            %                             if length(sigsen) == 1 %%% select the top n tagging sensors
                            %                                break
                            %                             end
                        end
                    end
                end
                %%% finally get the sig occipital tagging sensors
                cohcoh = nanmean(coh.cohspctrm(sigchanid,:,:),1);
                eval(['TagCoh.' EpochType{mmm} '_Coh(:,:,1,ccc) = squeeze(cohcoh);']);
            end
        end
        save([PPath.FigPath DataSets{ddd} '_' num2str(sss)],'TagCoh','-v7.3');
        
        %     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%                3. Plotting figures                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     %%%================= plotting coh and pow for single subject
        %     %%%%%% coh
        %     figtitle = [DataSets{ddd} '_' EpochType{eee} '_' RunCond '_s' num2str(sss) '_coh'];
        %     if isempty(sigsens)
        %         figtitle = ['NoSig_' DataSets{ddd} '_' EpochType{eee} '_' RunCond '_s' num2str(sss) '_coh'];
        %     end
        %     h = figure('Name',figtitle,'color',[1 1 1]);
        %     for mmm = 1:length(EpochType)
        %         %%%% coh plot
        %         eval(['tmpdata = nanmean(TagCoh.' EpochType{mmm} '_Coh(tagfreq_id_30,:,1,:),3);']);%%% freq*time*sub*cond
        %         tmpcc = [1 2 3];
        %         for cc = 1:length(conds)+1
        %             if ismember(cc,tmpcc)
        %                 condid = find(tmpcc==cc);
        %                 eval(['rt = TagCoh.' EpochType{mmm} '_RT(:,condid)./1000;']);
        %                 rt = nanmean(rt);
        %                 loc = tmpcc(condid)+4*(mmm-1);
        %                 plotdata = squeeze(tmpdata(:,:,:,condid));
        %                 subtitle = [EpochType{mmm} 'et  ' TagCoh.conds{condid}];
        %                 subplot(3,4,loc);
        %                 pcolor(TagCoh.time, TagCoh.freq(tagfreq_id_30), plotdata); colorbar;
        %                 %                 caxis([0 .4]);
        %             else
        %                 eval(['rt = TagCoh.' EpochType{mmm} '_RT(:,1)./1000;']); %% for the difference figure, use the RT from 'all' cond
        %                 rt = nanmean(rt);
        %                 loc = cc+4*(mmm-1);
        %                 if cc == 4
        %                     idx1 = 2; idx2 = 3;
        %                 end
        %                 plotdata = squeeze((tmpdata(:,:,:,idx1)-tmpdata(:,:,:,idx2))./(tmpdata(:,:,:,idx1)+tmpdata(:,:,:,idx2)));
        %                 subtitle = [EpochType{mmm} 'et Coh-contrast'];
        %                 subplot(3,4,loc);
        %                 pcolor(TagCoh.time, TagCoh.freq(tagfreq_id_30), plotdata); colorbar;
        %                 %                 caxis([-.07 .07]);
        %             end
        %             colormap jet; shading interp; % colorbar;
        %             hold on;
        %             plot([-0.5 0.5],[60 60],'-.k','LineWidth',2)
        %             plot([-rt -rt],[1 100],'-.k','LineWidth',2)
        %             plot([0 0],[1 100],'-.k','LineWidth',2)
        %             plot([rt rt],[1 100],'-.k','LineWidth',2)
        %             set(gca,'XTick',-0.4:0.2:0.4);
        %             set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'},'FontWeight','bold','FontSize',10)
        %             title(subtitle,'FontWeight','bold','FontSize',14);
        %             %             if mmm == 1 && (cc == 1 || cc ==4 )
        %             %                 colorbar('FontWeight','bold','FontSize',10)
        %             %             end
        %         end
        %     end
        %     %     saveas(h,[PPath.FigPath figtitle]);
        %     ScSz = [0 0 1300 800];
        %     set(gcf,'Position',ScSz);
        %     print('-dpng', '-r300', [PPath.FigPath figtitle])
        %     close all
        %
    end
    toc
    
else  %%%server == 0
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run locally get the whole structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% basic parameters
    %%% combine data together into one big structure
    load([PPath.FigPath DataSets{ddd} '_3']);
    %%% since the first sub is sub_3, but sub_3's TagCoh has no subject
    %%% dimension, which is index 1 actually, so need to remove the data,
    %%% otherwise the sub_3 data will show up in the sub_1 slot
    TagCoh_all = TagCoh;
    nsub = length(TagCoh.subs);
    TagCoh_all.SigTagSen_Occip(1:nsub,1) = repmat({'No'},nsub,1);
    TagCoh_all.SigTagSen_Occip(1:nsub,2) = repmat({0},nsub,1);
    for ps = sss_targ'
        load([PPath.FigPath DataSets{ddd} '_' num2str(ps)]);
        for fd = 1:length(EpochType)
            if ps == 1
                TagCoh_all.Stat(:,ps) = [];
                eval(['TagCoh_all.' EpochType{fd} '_TrlNum_raw(ps,:)= [];']);
                eval(['TagCoh_all.' EpochType{fd} '_RT(ps,:)= [];']);
                eval(['TagCoh_all.' EpochType{fd} '_TrlInfo(ps,:)= [];']);
                eval(['TagCoh_all.' EpochType{fd} '_TrlInfo(ps,:)= [];']);
                TagCoh_all.SigTagSen_All{ps,1} = []; %%% all sig tagging occipital sensors
                eval(['TagCoh_all.' EpochType{fd} '_Coh(:,:,ps,:)= [];']);
            end
            TagCoh_all.Stat(:,ps) = TagCoh.Stat(:,1);
            eval(['TagCoh_all.' EpochType{fd} '_TrlNum_raw(ps,:)= TagCoh.' EpochType{fd} '_TrlNum_raw(1,:);']);
            eval(['TagCoh_all.' EpochType{fd} '_RT(ps,:)= TagCoh.' EpochType{fd} '_RT(1,:);']);
            eval(['TagCoh_all.' EpochType{fd} '_TrlInfo(ps,:)= TagCoh.' EpochType{fd} '_TrlInfo(1,:);']);
            TagCoh_all.SigTagSen_All{ps,1} = TagCoh.SigTagSen_All{1,1}; %%% all sig tagging occipital sensors
            TagCoh_all.SigTagSen_Occip(ps,[1 2]) = TagCoh.SigTagSen_Occip(1,[1 2]);
            eval(['TagCoh_all.' EpochType{fd} '_Coh(:,:,ps,:)= TagCoh.' EpochType{fd} '_Coh(:,:,1,:);']);
        end
    end
    clear TagCoh
    TagCoh = TagCoh_all;
    clear TagCoh_all
    TagCoh.SenSelectP = 0.01;
    lesstrl_sub = find(sum(TagCoh.PreTarg_TrlNum_raw(:,[2 3])<60,2)>0); %%% including non-sig subs
    nosigsub = find(strcmp(TagCoh.SigTagSen_Occip,'No'));
    TagCoh.LessTrlSubID = lesstrl_sub;
    TagCoh.SigSubID = transpose(setdiff(1:length(TagCoh.subs),[nosigsub; lesstrl_sub]));
    TagCoh.MeanPretargTrlNum4SigSub = nanmean(TagCoh.PreTarg_TrlNum_raw(TagCoh.SigSubID,:));
    save([PPath.FigPath DataSets{ddd}],'TagCoh','-v7.3');
    % delete([PPath.FigPath DataSets{ddd}  '_*']);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. group plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plotting-coh
    sigsubid = TagCoh.SigSubID;
    figtitle = [DataSets{ddd} '_' EpochType{eee} '_' RunCond '_Coh'];
    h = figure('Name',figtitle,'color',[1 1 1]);
    ScSz = [1 1 900 800];%get(groot, 'Screensize' );
    set(h,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
    pos_grid =  [ScSz(1)+10 ScSz(2)+10  ScSz(3)-10  ScSz(4)-10];
    nrow = length(conds)+1;
    ncol = length(EpochType);
    gridPos = DivideScreen(nrow,ncol,ScSz,50,70,pos_grid);
    for mmm = 1:length(EpochType)
        %%%% coh plot
        eval(['tmpall = TagCoh.' EpochType{mmm} '_Coh(:,:,sigsubid,:);']);%%% freq*time*sub*cond
        tmpdata = squeeze(nanmean(tmpall,3));%%% freq*time*sub*cond
        tmpcc = [1 2 3];
        for cc = 1:nrow
            if ismember(cc,tmpcc)
                condid = find(tmpcc==cc);
                eval(['rt = TagCoh.' EpochType{mmm} '_RT(sigsubid,condid)./1000;']);
                rt = nanmean(rt);
                loc = mmm + ncol*(tmpcc(condid)-1);
                plotdata = squeeze(tmpdata(:,:,condid));
                subtitle = [EpochType{mmm} 'et  ' TagCoh.conds{condid}];
                h = axes('Position',gridPos{loc});
                pcolor(TagCoh.time, TagCoh.freq, plotdata);axcopy;
                caxis([0 0.05]);
                xlim([-0.4 0.4])
            else
                loc = mmm + ncol*(cc-1);
                if cc == 4
                    idx1 = 2; idx2 = 3;
                else
                    idx1 = 4; idx2 = 5;
                end
                eval(['rt = TagCoh.' EpochType{mmm} '_RT(sigsubid,[idx1 idx2])./1000;']); %% for the difference figure, use the RT from 'all' cond
                rt = nanmean(rt,1);
                rt = min(rt);
%                 plotdata = squeeze(nanmean((tmpall(:,:,:,idx1)-tmpall(:,:,:,idx2))./(tmpall(:,:,:,idx1)+tmpall(:,:,:,idx2)),3));
                plotdata = squeeze(nanmean(tmpall(:,:,:,idx1)-tmpall(:,:,:,idx2),3));
                subtitle = [EpochType{mmm} 'et low-high'];
                h = axes('Position',gridPos{loc});
                pcolor(TagCoh.time, TagCoh.freq, plotdata);axcopy;
%                 caxis([-0.006 0.006]);
                xlim([-0.4 0.4])
            end
            colormap jet; shading interp;
            hold on;
            plot([-0.5 0.5],[60 60],'-.k','LineWidth',2)
            plot([-rt -rt],[1 100],'-.k','LineWidth',2)
            plot([0 0],[1 100],'-.k','LineWidth',2)
            plot([rt rt],[1 100],'-.k','LineWidth',2)
            set(gca,'XTick',-0.4:0.2:0.4);
            set(gca,'XTickLabel',{'-0.4','-0.2',['Fix' RunCond(4:end)],'0.2','0.4'},'FontWeight','bold','FontSize',12)
            title(subtitle,'FontWeight','bold','FontSize',16);
            if mmm == 3 && (cc == 1 || cc ==4 )
                colorbar('FontWeight','bold','FontSize',10)
            end
        end
    end
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'bmp');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Get total topograph of RFT %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% get the occipital sensors
    load([rootdir filesep 'Analyse_data' filesep '20190801_b32f' filesep 'Targ60_TW1000' filesep 'epoch_BL_Cross.mat']) %%% just load random data that has all the labels
    cfg = [];
    erf = ft_timelockanalysis(cfg, epoch_BL_Cross);
    
    %%% get labels
    figtitle = [DataSets{ddd} '-' EpochType{eee} '-' RunCond '-RFT-sens'];
    SigSen = TagCoh.SigTagSen_Occip(sigsubid,1);
    tmp = zeros(size(erf.label));
    for bbb = 1:length(SigSen)
        for kkk = 1:length(SigSen{bbb})
            idx = find(strcmp(erf.label,[SigSen{bbb}{kkk}(1:end-1) '1']));
            tmp(idx) = tmp(idx)+1;
        end
    end
    erf.avg = repmat(tmp,1,size(erf.avg,2));
    h = figure('color', [1 1 1],'name',figtitle);
    cfg =[];
    cfg.layout = 'neuromag306mag.lay'; %'neuromag306cmb.lay';
    cfg.comment = ' ';
    %     cfg.gridscale = 300;
    ft_topoplotER(cfg,erf);
    colormap jet; colorbar('FontWeight','bold','FontSize',10);
    caxis([-max(tmp) max(tmp)])
    title('RFT reponse sensors','FontWeight','bold','FontSize',16);
    text(0.5,0.5,'No. of sensors','FontWeight','bold','FontSize',14,'Rotation',270)
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'bmp');
    % close all
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Just do the simple ttest to find the effect %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq60 = dsearchn(TagCoh.freq',60);
    sigsubid = TagCoh.SigSubID;
    nsub = length(sigsubid);
    fig_sufix = {''}; %%{'','TFR'};
    figtitle = [DataSets{ddd}  '_' RunCond '_Ttest_Freq'];
    h = figure('Name',figtitle,'color',[1 1 1]);
    for fff = 1 %:length(fig_sufix)
        eval(['tmptime = transpose(TagCoh.time' fig_sufix{fff} ');']);
        t0 = dsearchn(tmptime,0);
        for mm = 1:length(EpochType)
            eval(['rt = TagCoh.' EpochType{mm} '_RT;']);
            eval(['rt = TagCoh.' EpochType{mm} '_TrlInfo;']);
            if isempty(fig_sufix{fff})
                eval(['cohdata = TagCoh.' EpochType{mm} '_Coh;']);
            else
                eval(['cohdata = TagCoh.' EpochType{mm} '_PowRaw;']);
            end
            tdata = zeros(nsub,2);
            trl_minrt = zeros(nsub,1);
            for ss = 1:nsub
                sid = sigsubid(ss);
                % %                 %%% get min rt of 2 averaged rt in 2 conditions
                % %                 t1 = dsearchn(tmptime,rt(sid,2)/1000);
                % %                 t2 = dsearchn(tmptime,rt(sid,3)/1000);
                % %                 ttt = min([t1 t2]);
                %%% get min rt of all the trials in 2 conditions
                trl_minrt(ss,1) = min([rt{sid,2}(:,8); rt{sid,3}(:,8)])/1000;
                if strfind(RunCond,'On') %% epoch aligned with fixaiton_on, select timewindow from zero to this fixaiton_duration
                    ttt = dsearchn(tmptime,trl_minrt(ss,1));
                    tdata(ss,:) = squeeze(nanmean(cohdata(freq60,t0:ttt,sid,[2 3]),2));
                elseif strfind(RunCond,'Off') %% epoch aligned with fixaiton_off, select timewindow from this fixaiton_duration to zero
                    ttt = dsearchn(tmptime,-trl_minrt(ss,1));
%                     ttt = dsearchn(tmptime,-0.2);
                    tdata(ss,:) = squeeze(nanmean(cohdata(freq60,ttt:t0,sid,[2 3]),2)); 
                end
            end
            [~,p,~,stats] = ttest(tdata(:,1), tdata(:,2));
            stats.p = p; disp([EpochType{mm} '_Ttest' fig_sufix{fff} '= ' num2str(p)])
            stats.data4test = tdata;
            stats.trl_minrt = trl_minrt;
            eval(['TagCoh.' EpochType{mm} '_Ttest' fig_sufix{fff} '= stats;'])
            %%% plot
            subplot(length(fig_sufix),3,mm+(fff-1)*3)
            x1 = 1.*ones(nsub,1); x2 = 2.*ones(nsub,1);
            y1 = tdata(:,1); y2 = tdata(:,2);
            hold on;
            scatter(x1,y1,50,[.8 .8 .8]);
            scatter(x2,y2,50,[.8 .8 .8]);
            plot([x1 x2]',[y1 y2]','Color',[.8 .8 .8])
            %%% plot averaged data
            scatter(x1(1),nanmean(y1),50,[1 0 0],'filled');
            scatter(x2(1),nanmean(y2),50,[1 0 0],'filled');
            plot([x1(1) x2(1)]',[nanmean(y1) nanmean(y2)]','Color',[1 0 0],'LineWidth',2)
            if isempty(fig_sufix{fff})
                title([EpochType{mm} ' Coh p=' num2str(round(1000*p)/1000)],'FontSize',10);
            else
                title([EpochType{mm} ' Pow p=' num2str(round(1000*p)/1000)],'FontSize',10);
            end
            set(gca,'xticklabel',{'Low','','High'})
        end
    end
    saveas(h,[PPath.FigPath figtitle]);
    
    %% %%%======= violin plot of the word-freq effect on coherence
    colmat = [0 114 189;217 83 25]./255;
    figtitle = [DataSets{ddd}  '_' RunCond '_Ttest_Freq_violin'];
    nsigsub = size(TagCoh.PreTarg_Ttest.data4test,1);
    group = [cellstr(repmat('Low',nsigsub,1)); cellstr(repmat('High',nsigsub,1))];
    grouporder={'Low','High'};
    figure('Name',figtitle,'color',[1 1 1],'Position',[1 1 1000 300]);
    for eee = 1:length(EpochType)
        eval(['vdata = [TagCoh.' EpochType{eee} '_Ttest.data4test(:,1); TagCoh.' EpochType{eee} '_Ttest.data4test(:,2)];']);
        h = subplot(1,length(EpochType),eee);
        vp = violinplot(vdata, group,'GroupOrder',grouporder);
        vp(1).ViolinColor = colmat(1,:);
        vp(2).ViolinColor = colmat(2,:);
        ylabel('Coherence (r)','FontSize',16,'FontWeight','bold');
        set(gca,'FontSize',12,'FontWeight','bold');
        set(gca,'box','off','LineWidth',2)
        yyy = get(h,'ylim');
        eval(['pvalue = TagCoh.' EpochType{eee} '_Ttest.p;']);
        title([EpochType{eee} 'et' 'p = ' num2str(round(1000*pvalue)/1000)],'FontSize',16,'FontWeight','bold');
        %%% plot the line linking each subject
        hold on;
        x1 = vp(1,1).ScatterPlot.XData;
        y1 = vp(1,1).ScatterPlot.YData;
        x2 = vp(1,2).ScatterPlot.XData;
        y2 = vp(1,2).ScatterPlot.YData;
        plot([x1; x2],[y1; y2],'Color',[.7 .7 .7])
    end
    xlabel('Word Frequency for Target','FontSize',16,'FontWeight','bold');
    saveas(h,[PPath.FigPath figtitle]);
    % %     saveas(h,[PPath.FigPath figtitle],'bmp');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% coherence latency statistics ---- JakeNife %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colmat = [0 114 189;217 83 25]./255;
    condtype = {'WrdFreq'};
    legend_all = [{[]},{'low'},{'high'}];
    HM_tim = zeros(length(sigsubid),length(EpochType),length(conds)); %%% sub*epoch*condition
    PK_tim = HM_tim; %%% sub*epoch*condition
    HM_tim_grd = zeros(length(EpochType),length(conds)); %%% epoch*condition
    PK_tim_grd = HM_tim_grd;
    %%% latency and magnitude of 'center of mass'
    CM_tim = HM_tim;
    CM_tim_grd = HM_tim_grd;
    sigsubid = TagCoh.SigSubID;
    tagfreqid = 11;
    %%% get the time range for important timepoints
    timstp = 0.001;
    min_tp = 50*timstp; %%% 50ms
    %%% plot RFT curve for all conditions
    for i = 1:length(EpochType)
        eval(['tempdata = TagCoh.' EpochType{i} '_Coh;']);
        tempdata = squeeze(nanmean(tempdata(tagfreqid,:,sigsubid,:),1)); %%% time*sub*cond
        eval(['TagCoh.' EpochType{i} '_Coh_TimSubCond = tempdata;']);
    end
    
    epoch_type = {'PreTarg','Targ'};  %;'PosTarg'};
    for ee = 1:length(epoch_type)
        eval(['cohdata_all = TagCoh.' epoch_type{ee} '_Coh_TimSubCond(:,:,[2 3]);']); %% time*sub*cond
        eval(['RT_all = TagCoh.' epoch_type{ee} '_RT(sigsubid,[2 3])./1000;']); %% sub*cond
        CM_tim = []; %% center of mass
        HM_tim = []; %% when reach the half maximum amplitude (distant amplitude between max and min coh)
        
        % %         %%% parameters for plotting --- each iteration of the jack-knife
        % %         figtitle = [DataSets{ddd} '_' epoch_type{ee} '_Latency_JkNf_allsub'];
        % %         h = figure('Name',figtitle,'color',[1 1 1]);
        % %         ScSz = [1 1 1400 600];%get(groot, 'Screensize' );
        % %         set(h,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
        % %         pos_grid =  [ScSz(1)+5 ScSz(2)+5  ScSz(3)-5  ScSz(4)-5];
        % %         gridPos = DivideScreen(6,ceil(length(sigsubid)/6),ScSz,35,90,pos_grid);
        
        %%% doing jackknife loop
        allsigsub = 1:length(sigsubid);
        exclu_sub = [allsigsub 0]; %%[jackknife all]
        for sid = 1:length(exclu_sub)
            subid_jn = setdiff(allsigsub,exclu_sub(sid));
            %%% data after jackknife
            cohdata = cohdata_all(:,subid_jn,:);
            RT = RT_all(subid_jn,:);
            meanRFT = squeeze(nanmean(cohdata,2)); %%% time*cond
            seRFT = squeeze(nanstd(cohdata,0,2)./sqrt(size(cohdata,2)));
            meanRT = mean(RT);
            %%% get the time window
            rt = min(meanRT);
            etp = find(TagCoh.time < rt); %%% rt window aligned with zero--saccadeonset
            zero_tp = nearest(TagCoh.time,0);
            tprange = zero_tp:etp(end);
            timrange = TagCoh.time(tprange);
            
            for cc = 1:2 %%[low high]
                %%%=== center of mass
                tmptp = meanRFT(tprange,cc); %%% data in the timewindow: time*sub*cond
                timevct = [1:length(timrange)]';
                cm = round(sum(timevct.*tmptp)/sum(tmptp));
                CM_tim(sid,cc) = timrange(cm);
                %%% percentage-amplitude
                %             percent = 1;
                %             tmptp = tmptp-min(tmptp);
                %             halfcoh_id = nearest(tmptp,percent*max(tmptp));
                %             HM_tim(sid,cc,ee) = timrange(halfcoh_id);
                %%%=== half amplitude latency
                [mintmp, minidx] = min(tmptp);
                %%% make sure that the maxtmp is later than the mintmp
                tmptp = tmptp(minidx:end);
                [maxtmp, maxidx] = max(tmptp);
                tmpidx = find(tmptp>((maxtmp-mintmp)/2+mintmp));
                if ~isempty(tmpidx)
                    HM_tim(sid,cc) = timrange(tmpidx(1)+minidx);
                else
                    HM_tim(sid,cc) = nan;
                end
            end
            
            % %             %%% plot to check
            % %             hhh = axes('Position',gridPos{sid});
            % %             rt = 0.4; %%% just for the plot
            % %             etp = find(TagCoh.time < rt); %%% rt window aligned with zero--saccadeonset
            % %             zero_tp = nearest(TagCoh.time,0);
            % %             tprange = zero_tp:etp(end);
            % %             timrange = TagCoh.time(tprange);
            % %             hold on;
            % %             a = shadedErrorBar(timrange,meanRFT(tprange,1),seRFT(tprange,1),{'color',colmat(1,:)},0.5);
            % %             b = shadedErrorBar(timrange,meanRFT(tprange,2),seRFT(tprange,2),{'color',colmat(2,:)},0.5);
            % %             timeliney = get(hhh,'ylim');
            % %             plot([HM_tim(sid,1,ee) HM_tim(sid,1,ee)],timeliney,'color',colmat(1,:)); %%% low
            % %             plot([HM_tim(sid,2,ee) HM_tim(sid,2,ee)],timeliney,'color',colmat(2,:)); %%% high
            % %             axcopy;
            % %             legendflex([a.mainLine,b.mainLine],[{'Low'},{'High'}],'anchor', {'ne','ne'}, 'buffer', [5 -5],'FontWeight','bold','FontSize',10,'xscale',0.8,'box', 'off');
            % %             if sid == length(exclu_sub)
            % %                 title([epoch_type{ee} '-HM-Group'],'FontWeight','bold','FontSize',10)
            % %             else
            % %                 title([epoch_type{ee} '-HM-NoSub-' num2str(sid)],'FontWeight','bold','FontSize',10)
            % %             end
        end
        % %         saveas(h,[PPath.FigPath figtitle]);
        eval(['JkNf.' epoch_type{ee} '_HM_tim = HM_tim ;']);
        eval(['JkNf.' epoch_type{ee} '_CM_tim = CM_tim ;']);
    end
    
    %%% using jacknife to do the statistics
    tp_tpye = {'CM_tim';'HM_tim'};
    N = length(sigsubid);
    for tp = 1:length(tp_tpye)
        for ee = 1:length(epoch_type)
            tj = [];
            eval(['tmptim = JkNf.' epoch_type{ee} '_' tp_tpye{tp} ';']);
            Di = tmptim(1:end-1,1)-tmptim(1:end-1,2);
            J = nansum(Di)/N;
            SD = sqrt([(N-1)/N]*[nansum((Di-J).^2)]);
            D = tmptim(end,1) - tmptim(end,2);
            tj = D/SD;
            eval(['JkNf.'  epoch_type{ee} '_' tp_tpye{tp} '_tvalue = tj;']);
        end
    end
    TagCoh.JkNf = JkNf;
    save([PPath.FigPath DataSets{ddd} '.mat'], 'TagCoh')
    
    %% === plot the group level curve with sig, only for HM (half-maximum)
    curve_xmax = 0.4; %%% the x-axis range
    etp = find(TagCoh.time < curve_xmax); %%% rt window aligned with zero--saccadeonset
    zero_tp = nearest(TagCoh.time,0);
    tprange = zero_tp:etp(end);
    timrange = TagCoh.time(tprange);
    cohdata_all = TagCoh.PreTarg_Coh_TimSubCond(tprange,:,[2 3]); %% time*sub*cond
    meanRFT = squeeze(mean(cohdata_all,2));
    seRFT = squeeze(nanstd(cohdata_all,0,2)./sqrt(length(sigsubid)));
    HM = TagCoh.JkNf.HM_tim(end,:,1);
    colmat = [0 114 189;217 83 25]./255;
    %%%% figure
    figtitle = 'PreTarget coherence onset latency';
    h = figure('Name',figtitle,'color',[1 1 1]);
    a = shadedErrorBar(timrange,meanRFT(:,1),seRFT(:,1),{'color',colmat(1,:)},0.8);hold on;
    b = shadedErrorBar(timrange,meanRFT(:,2),seRFT(:,2),{'color',colmat(2,:)},0.9);
    a.mainLine.LineWidth = 2;
    b.mainLine.LineWidth = 2;
    legendflex([a.mainLine,b.mainLine],{'low';'high'},'anchor', {'ne','ne'}, 'buffer', [0 0],'FontWeight','bold','Fontsize',12,'xscale',1,'box','off');
    % % %     %%%% no sig for 2 datasets cmb
    % % %     text(0.3,0.04,'n.s.','FontSize',20);
    % %     %%%%% sig for JEP60 dataset
    % %     text(0.3,0.05,'p < 0.01','FontWeight','bold','FontSize',12);
    text(0.3,0.04,['t-value is ' num2str(round(1000*TagCoh.JkNf.HM_tim_tvalue)/1000)],'FontWeight','bold','FontSize',12);
    set(gca,'FontSize',12,'FontWeight','bold');
    set(gca,'box','off','LineWidth',2)
    set(gca,'XTick',0:0.05:0.4);
    set(gca,'XTickLabel',{'FixOn','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4'},'FontWeight','bold','FontSize',12)
    set(gca,'YTick',0.01:0.01:0.06);
    title(figtitle,'FontWeight','bold','FontSize',16)
    xlabel('Time (ms)','FontWeight','bold','FontSize',16)
    ylabel('PreTarget coherence (r)','FontWeight','bold','FontSize',16)
    %%% plot the latency lines
    hold on;
    timeliney = get(gca,'ylim');
    timelinex = repmat(HM(1),1,2);
    plot(timelinex', timeliney','--','color',colmat(1,:),'LineWidth',2);
    timelinex = repmat(HM(2),1,2);
    plot(timelinex', timeliney','--','color',colmat(2,:),'LineWidth',2);
    saveas(h,[PPath.FigPath figtitle]);
    % %     saveas(h,[PPath.FigPath figtitle],'tiff');
    
    
    % %     %% %% plotting individual curve and mark the important time points in pre-target and target periods for frequency condition
    % %     type_loop = {'PreTarg';'Targ'};
    % %     colmat = [0 114 189;217 83 25]./255;
    % %     for fff = 1%:length(type_loop)
    % %         eval(['tempdata = TagCoh.' type_loop{fff} '_Coh_TimSubCond;']);
    % %         eval(['RT_all = TagCoh.' type_loop{fff} '_RT(sigsubid,:)./1000;']); %% sub*cond
    % %         meanRFT = squeeze(nanmean(tempdata,2)); %%% time*cond
    % %         seRFT = squeeze(nanstd(tempdata,0,2)./sqrt(length(sigsubid)));  %%% time*cond
    % %         %%% plotting
    % %         figtitle = [DataSets{ddd} '_' type_loop{fff} '_HM_subcurve'];
    % %         h = figure('Name',figtitle,'color',[1 1 1]);
    % %         ScSz = [1 1 1300 600];%get(groot, 'Screensize' );
    % %         set(h,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
    % %         pos_grid =  [ScSz(1)+20 ScSz(2)+5  ScSz(3)-5  ScSz(4)-5];
    % %         gridPos = DivideScreen(4,ceil(length(sigsubid)/4),ScSz,50,110,pos_grid);
    % %         for k1 = 1:length(sigsubid)
    % %             loc = k1;
    % %             hhh = axes('Position',gridPos{loc});
    % %             for k2 = 2:3 %% low, high
    % %                 rt = min(RT_all(k1,[2 3]));
    % %                 etp = find(TagCoh.time < rt);
    % %                 tprange = find(TagCoh.time==min_tp):etp(end);
    % %                 tmptp = tempdata(tprange,k1,k2); %%% time*sub*cond
    % %                 timrange = TagCoh.time(tprange);
    % %                 hold on;
    % %                 plot(timrange,tmptp,'LineWidth',2);
    % %                 title(['Sub-' num2str(k1) '-HM'],'FontWeight','bold','FontSize',10)
    % %                 set(gca,'FontSize',8,'FontWeight','bold');
    % %             end
    % %             if k1 == 1
    % %                 xlabel('Time (ms)','FontWeight','bold','FontSize',10)
    % %                 ylabel('PreTarget coh (r)','FontWeight','bold','FontSize',8)
    % %             end
    % %             axcopy;
    % %         end
    % % % % %         %%% plot group curve
    % % % % %         loc = k1+1;
    % % % % %         hhh = axes('Position',gridPos{loc});
    % % % % %         tprange = 501:900;
    % % % % %         for k2 = 2:3 %% low, high
    % % % % %             timrange = TagCoh.time(tprange);
    % % % % %             hold on;
    % % % % %             if k2 == 2
    % % % % %                 a = shadedErrorBar(timrange,meanRFT(tprange,k2),seRFT(tprange,k2),{'color',colmat(k2-1,:)},0.5);
    % % % % %             else
    % % % % %                 b = shadedErrorBar(timrange,meanRFT(tprange,k2),seRFT(tprange,k2),{'color',colmat(k2-1,:)},0.5);
    % % % % %             end
    % % % % %         end
    % % % % %         %%%% plot the important time point!
    % % % % %         timeliney = get(hhh,'ylim');
    % % % % %         plot([HM(1) HM(1)],timeliney,'color',colmat(1,:)); %%% low
    % % % % %         plot([HM(2) HM(2)],timeliney,'color',colmat(2,:)); %%% high
    % % % % %         axcopy;
    % % % % %         legendflex([a.mainLine,b.mainLine],[{'Low'},{'High'}],'anchor', {'ne','ne'}, 'buffer', [5 -5],'FontWeight','bold','FontSize',8,'xscale',0.8,'box', 'off');
    % % % % %         title([type_loop{fff} '-HM-Group'],'FontWeight','bold','FontSize',10)
    % % % % %         xlabel('Time (ms)','FontWeight','bold','FontSize',10)
    % % % % %         ylabel('PreTarget coh (r)','FontWeight','bold','FontSize',8)
    % % % % %         set(gca,'FontSize',8,'FontWeight','bold');
    % %         saveas(h,[PPath.FigPath figtitle]);
    % %         saveas(h,[PPath.FigPath figtitle],'tiff');
    % %     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%=== correlation between RFT and behavioral data =====%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %     %%% get the combined beha data
    % % %     sigsubid = TagCoh.SigSubID;
    % % %     sigsubname = TagCoh.subs(sigsubid);
    % % %     load([rootdir 'Analyse_data' filesep 'SubInfo.mat']);
    % % %     Targ60_BD = load([rootdir filesep 'Results' filesep 'Behavioral' filesep 'Targ60_BehaData']);
    % % %     JEP60_BD = load([rootdir filesep 'Results' filesep 'Behavioral' filesep 'JEP60_BehaData']);
    % % %     %%%+++ get the subs that have 2 datasets
    % % %     subtarg = SubInfo.subjects.Targ60;
    % % %     subjep = SubInfo.subjects.JEP60;
    % % %     BehaData = [];
    % % %     loop_1 = {'FirstFix','Gaze','TotalGaze','PupilSize'};
    % % %     loop_2 = {'Pre','Tag','Pos'};
    % % %     for tmpsub = 1:length(subjep)  %%% sub id in JEP60
    % % %         stmp = find(strcmp(subjep{tmpsub}(1:end-4),subtarg)); %%% sub id in Targ60
    % % %         if stmp
    % % %             for p1 = 1:length(loop_1)
    % % %                 for p2 = 1:length(loop_2)
    % % %                     eval(['BehaData.' loop_1{p1} '.' loop_2{p2} '(stmp,:) = mean([Targ60_BD.BehaData.' loop_1{p1} '.' loop_2{p2} '(stmp,:); JEP60_BD.BehaData.' loop_1{p1} '.' loop_2{p2} '(tmpsub,:)],1);']);
    % % %                 end
    % % %             end
    % % %         end
    % % %     end
    % % %     save([rootdir filesep 'Results' filesep 'Behavioral' filesep 'CmbDatasets_BehaData'],'BehaData');
    
    % %     %%% load beha data
    % %     load([rootdir filesep 'Results' filesep 'Behavioral' filesep 'Targ60_BehaData']);
    % %     load([rootdir filesep 'Results' filesep 'Behavioral' filesep 'JEP60_BehaData']);
    load([rootdir filesep 'Results' filesep 'Behavioral' filesep 'CmbDatasets_BehaData']);
    
AllFix = nanmean(BehaData.FirstFix.Tag(sigsubid,:),2);
DiffFix_Targ = BehaData.FirstFix.Tag(sigsubid,:);
DiffFix_2 = DiffFix_Targ(:,1) - DiffFix_Targ(:,2); %% diff fixation for target: low-high
figtitle = [DataSets{ddd} '-Corr'];
h = figure('Name',figtitle,'color',[1 1 1]);
preTP = TagCoh.PreTarg_Ttest.data4test; %%% coherence of pretarget
DiffTP = preTP(:,1) - preTP(:,2); %%% low-high
%%%% PreTarg-DiffRFT & Targ Diff Fixation Duration
subplot(1,3,1)
scatter(DiffTP,DiffFix_2,50,'r','filled')
ylabel('Target First Fixation low-high (ms)','FontWeight','bold','FontSize',16)
xlabel('PreTarget coherence low-high (r)','FontWeight','bold','FontSize',16)
set(gca,'FontSize',12,'FontWeight','bold');
set(gca,'box','off','LineWidth',2)
[coef, pval] = corr([DiffTP  DiffFix_2],'type','Spearman');
title(['r = ' num2str(round(1000*coef(1,2))/1000) ', p = ' num2str(round(1000*pval(1,2))/1000)], 'FontWeight','bold','FontSize',16);
xlim([-0.02 0.03])
ylim([-30 70])

%%% PreTarg-DiffRFT &  Targ Fixation Duration
subplot(1,3,2)
scatter(DiffTP,AllFix,50,'r','filled')
ylabel('Target First Fixation (ms)','FontWeight','bold','FontSize',16)
xlabel('PreTarget coherence low-high (r)','FontWeight','bold','FontSize',16)
set(gca,'FontSize',12,'FontWeight','bold');
set(gca,'box','off','LineWidth',2)
[coef, pval] = corr([DiffTP AllFix],'type','Spearman');
title(['r = ' num2str(round(1000*coef(1,2))/1000) ', p = ' num2str(round(1000*pval(1,2))/1000)], 'FontWeight','bold','FontSize',16);
xlim([-0.02 0.03])
ylim([100 450])

%%% PreTarg-DiffRFT & total RT
subplot(1,3,3)
sigsubid = TagCoh.SigSubID;
SentDur = BehaData.SentDur_perwrd(sigsubid);
scatter(DiffTP,SentDur,50,'r','filled')
ylabel('Averaged fxiation in sentence (ms)','FontWeight','bold','FontSize',12)
xlabel('PreTarget coherence low-high (r)','FontWeight','bold','FontSize',12)
set(gca,'FontSize',12,'FontWeight','bold');
set(gca,'box','off','LineWidth',1.5)
[coef, pval] = corr([DiffTP  SentDur],'type','Spearman');
title(['r = ' num2str(round(1000*coef(1,2))/1000) ', p = ' num2str(round(1000*pval(1,2))/1000)], 'FontWeight','bold','FontSize',16);
xlim([-0.02 0.03])
ylim([50 360])
saveas(h,[PPath.FigPath figtitle]);    
save([PPath.FigPath DataSets{ddd} '.mat'], 'TagCoh')


SentDur_nopre = BehaData.SentDur_perwrd_nopre(sigsubid);
[coef, pval] = corr([DiffTP  SentDur_nopre],'type','Spearman') %%[-0.3258 0.1046]


end

% % % %%%% unifying plot properties:
% % % title/big x or y label: 'FontSize',16,'FontWeight','bold'
% % % text: * -- 20, n.s.--12
% % % axis: set(gca,'FontSize',12,'FontWeight','bold');


% % %% %%%================== run 2 factors ANOVA for first fixation ====================%%%%%
% % load('Z:\Results\Coh_hilbert_01\minirt_trl\Cmb.mat')
% % X = [TagCoh.PreTarg_Ttest.data4test; TagCoh.Targ_Ttest.data4test];
% % [p,Table,stats] = anova2(X,size(TagCoh.PreTarg_Ttest.data4test,1));
% % TagCoh.Pre_Tag_ANOVA2 = Table;
% % save('Z:\Results\Coh_hilbert_01\minirt_trl\Cmb.mat','TagCoh')


%%% TagCoh.PreTarg_Ttest.trl_minrt = 0.0883 +- 0.0089 ms (mean +- std)
%%% TagCoh.Targ_Ttest.trl_minrt = 0.0871 +- 0.0094 ms (mean +- std)










