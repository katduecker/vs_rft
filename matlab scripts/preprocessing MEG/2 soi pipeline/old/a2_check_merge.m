%% VS + RFT
% PhD project 2
% [c] Katharina Duecker

% script used to further investigate inconsistencies identified by
% a_merge_edf_fif_mat


clear all; close all; clc; beep off
% define paths
pth = 'Z:\Visual Search RFT';
%mtpth = fullfile('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT','matlab scripts');
addpath('Z:\fieldtrip')
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
mkdir(trl_merge_pth)
ft_defaults;

% list subj
d = dir(trl_merge_pth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds


load(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))
% discard rejected subj
subjfolds = subjfolds(ismember(subjfolds,mergesubj(:,1)));

% double check and align merge documentation
for s = 1:length(subjfolds)
    if ~isempty(mergesubj{s,2})
        if exist(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
                'trl_overlap_meg_el_rsp.mat'))
            
            load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
                'trl_overlap_meg_el_rsp.mat'))
        else
            mergesubj{s,3} = 'not aligned';
            mergesubj{s,4} = 'not aligned';
            continue
        end
        
        [trials_ok, rt_ok, idx_rt_diff] = a1_check_merged_trls(rspinfo, meginfo);
        
        if ~trials_ok || ~rt_ok
            mergesubj{s,3} = 'trials and RT mismatch';
        else
            mergesubj{s,3} = 'adjusted';
            mergesubj{s,4} = 'adjusted';
            
        end
        
        ignore_diff = 1;
        save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
            'trl_overlap_meg_el_rsp.mat'),'rspinfo','elinfo','meginfo','idx_rt_diff','ignore_diff','-v7.3')
        
        save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_clean.mat'),'mergesubj')
        clear *info
    end
end

% check subjects RT: is the difference larger than 14 ms? (That's the
% trigger delay)
mergesubj(~ismember(mergesubj(:,1),subjfolds),:) = [];
idx_rt_diff_subj = ismember(mergesubj(:,3), 'trials and RT mismatch');
subj_rt_mismatch = subjfolds(idx_rt_diff_subj);
for s = 1:length(subj_rt_mismatch)
    load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
            'trl_overlap_meg_el_rsp.mat'))
        
        rsp_rt = [rspinfo.trl{:,3}]';
        rsp_rt(rspinfo.rsp_rejtrl) = [];
        diff_rt = meginfo.meg_rt-rsp_rt;
        
        % if the difference between MEG and stim computer is larger 14 ms,
        % discard
        diff_rt_large = diff_rt(abs(diff_rt)>0.014);
        
        % if all PTB RTs are > MEG RT
        if length(find(diff_rt_large>0)) == length(diff_rt_large)
            % you can ignore the difference and use the MEG RT
            ignore_diff = 1
            mergesubj{s,3} = 'adjusted, use MEG RT';
        else 
            idx_rt_diff = find(diff_rt>0.014);
            ignore_diff = 0
            mergesubj{s,3} = 'adjusted, use MEG RT, reject idx_rt_diff';
        end
        save(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},...
            'trl_overlap_meg_el_rsp.mat'),'rspinfo','elinfo','meginfo','idx_rt_diff','ignore_diff','-v7.3')
        
        
        save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_clean.mat'),'mergesubj')
end
