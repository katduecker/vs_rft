%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b1: RT & hit rate for alpha high vs low
% export .csv for statistical analysis in R

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Behavioural analyses
% a: performance per condition
% b. performance for alpha high vs low


%% Paths
clear all; close all; clc

split_ta_tp = 1;            % split for target present/absent?
toi_alpha_split = [-1 0];   % alpha toi

addpath('Z:\fieldtrip')
ft_defaults;

pth = 'Z:\Visual Search RFT';
addpath('Z:\Visual Search RFT\Violinplot-Matlab-master')
behavpth = fullfile(pth,'results','behavior');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','alpha high low');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;
toi_alpha_split = arrayfun(@num2str,toi_alpha_split.*1000,'UniformOutput',false);


if split_ta_tp

    condi = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},{'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};
    split_suf = '_ta_tp';


else
    condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
    split_suf = '';
    varnames = {'ung 16','gui 16', 'ung 32','gui 32'};

end
rt_alpha_high = zeros(length(subj),length(condi));
rt_alpha_low = rt_alpha_high;
% d_prime_alpha_high = rt_alpha_high;
% d_prime_alpha_low = rt_alpha_high;

hit_rate_high = rt_alpha_high;
hit_rate_low = rt_alpha_high;

for s = 1:length(subj)
    for c = 1:length(condi)

        % load results from coherence (don't re-run)
        load(fullfile(cohpth,subj{s},[strjoin(condi{c},'_'),'_fwidth',num2str(5),strjoin(toi_alpha_split,'_'),'_condi.mat']),'rt_high','rt_low','hit_high','hit_low')

        rt_alpha_high(s,c) = rt_high;
        rt_alpha_low(s,c) = rt_low;

        hit_rate_high(s,c) = hit_high;
        hit_rate_low(s,c) = hit_low;
 
    end
end

% code alpha high/low
high_low = zeros(length(subj)*length(condi)*2,1);
high_low(1:length(high_low)/2) = 1;

% condition 
condi_idx = arrayfun(@(x) repmat(x,length(subj),1),1:length(condi),'UniformOutput',false);
condi_idx = vertcat(condi_idx{:});
condi_idx = [condi_idx;condi_idx];


% code subj id
subj_id = repmat([1:length(subj)]',length(condi)*2,1);

RT_all = [reshape(rt_alpha_high,[],1);reshape(rt_alpha_low,[],1)];
hit_rate_all = [reshape(hit_rate_high,[],1);reshape(hit_rate_low,[],1)];


T = [subj_id,high_low,condi_idx,RT_all,hit_rate_all];
T = array2table(T,'VariableNames',{'id','high_low','condi','RT','hit rate'});

writetable(T,fullfile(behavpth,['RT_hitrate_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.csv']))

save(fullfile(behavpth,['RT_hitrate_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.mat']),'rt_alpha_high','rt_alpha_low','hit_rate_high','hit_rate_low')


rt_diff = mean(rt_alpha_high,2) - mean(rt_alpha_low,2);

avg_rt_diff = mean(rt_diff,2);


save(fullfile(behavpth,['RT_hitrate_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.mat']),'rt_alpha_high','rt_alpha_low','hit_rate_high','hit_rate_low','avg_rt_diff')

T = [mean(rt_alpha_high,2),mean(rt_alpha_low,2)];
T = array2table(T,'VariableNames',{'RT_high','RT_low'});

writetable(T,fullfile(behavpth,['RT_avg_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.csv']))

hit_diff = mean(hit_rate_high,2) - mean(hit_rate_low,2);

avg_hit_diff = mean(hit_diff,2);

save(fullfile(behavpth,['RT_hitrate_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.mat']),'avg_hit_diff','-append')

T = [mean(hit_rate_high,2),mean(hit_rate_low,2)]
T = array2table(T,'VariableNames',{'hit_high','hit_low'})

writetable(T,fullfile(behavpth,['hit_avg_alpha',strjoin(toi_alpha_split,'_'),split_suf,'.csv']))

