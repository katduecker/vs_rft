%% VS + RFT
% PhD project 2

% separate combined gradiometer labels

% [c] Katharina Duecker


pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
% matlab scripts path
mtpth = fullfile(pth,'matlab scripts','tfrs');
% load soi
load(fullfile(mtpth,'soi_tfr_subj.mat'))
% generate new cell in which grads are "pulled apart"

soigrad = cell(1,length(soicmb));
for i = 1:length(soicmb)
    cursoi = {};
    for l = 1:length(soicmb{i})
        cursoi = [cursoi, soicmb{i}{l}(1:7),['MEG',soicmb{i}{l}(9:end)]];
    end
    soigrad{i} = cursoi;
end

save(fullfile(mtpth,'soi_tfr_subj.mat'),'soicmb','soimag','soigrad')