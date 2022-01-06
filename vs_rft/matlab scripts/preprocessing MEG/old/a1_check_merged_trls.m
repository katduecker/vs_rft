%% VS + RFT
% PhD project 2

% check merged trials

% [c] Katharina Duecker

function [trials_ok, rt_ok,idx_rt_diff] = a1_check_merged_trls(rspinfo, meginfo)

rsptrl = [rspinfo.trl{:,1}]';
rsptrl(rspinfo.rsp_rejtrl,:) = [];
alltrl = meginfo.alltrl_list(:,1);
alltrl(meginfo.rejtrl_meg) = [];


if ~isequal(rsptrl,alltrl)
    warning('trials differ!')
    trials_ok = 0;
else
    trials_ok = 1;
end

% check RT
rsprt = [rspinfo.trl{:,3}]';
rsprt(rspinfo.rsp_rejtrl) = [];
% get rid of last 2 decimals

idx_rt_diff = find(abs(rsprt-meginfo.meg_rt)>0.014);

x = find(rsprt(idx_rt_diff) < meginfo.meg_rt(idx_rt_diff));

[rsprt(idx_rt_diff(x)),meginfo.meg_rt(idx_rt_diff(x))]
{rspinfo.trl{idx_rt_diff,2}}
{rspinfo.trl{idx_rt_diff(x),2}}

if ~isempty(idx_rt_diff)
    warning('RT''s differ!')
    rt_ok = 0;
else
    rt_ok = 1;
end

        
