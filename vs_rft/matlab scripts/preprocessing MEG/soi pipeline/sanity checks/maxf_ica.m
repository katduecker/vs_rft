%% Maxfilter vs ICA which subjects missing?

maxfpath = 'Z:\Visual Search RFT\results\meg\1 maxfilter';
mergepath = 'Z:\Visual Search RFT\results\meg\2 merged edf mat'; 
icapath = 'Z:\Visual Search RFT\results\meg\3 ICA';

d = dir(maxfpath);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));

d = dir(mergepath);
folds = {d.name};
subjfolds_merge = folds(strncmp(folds,'202',3));


d = dir(icapath);
files = {d.name};
subjfiles = files(strncmp(files,'202',3));
subjfiles = cellfun(@(x) x(1:13),subjfiles,'UniformOutput',false);

[diff_max_ica, pdiff] = setdiff(subjfolds,subjfiles);
[diff_merge_ica, pdiff_meica] = setdiff(subjfolds_merge,subjfiles);

pdiff_idx_max = find(ismember(subjfolds,subjfolds_merge(pdiff_meica)))

load(fullfile(mergepath,"docu_merge.mat"))
mergesubj(pdiff_idx_max,:)