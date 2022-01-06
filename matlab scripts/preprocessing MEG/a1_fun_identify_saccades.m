%% VS + RFT
% PhD project 2

% identify trials including saccades based on eyelink data

% [c] Katharina Duecker

function a1_fun_identify_saccades(s,pth,dtpth)
% server path
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');
addpath(fullfile(pth,'matlab scripts/','cbrewer'))

edfconvpath = fullfile(pth,'edf-converter-master/');
respth = fullfile(pth,'results/','eyelink/','1 identify saccades/');
trlpth = fullfile(pth,'results/','meg/','2 merged edf mat/');
pltpth = fullfile(respth,'plots');
mkdir(pltpth)


d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
load(fullfile(scriptpth,'idx_subjoi.mat'))

subj = subjfolds(usable_idx);

% screen settings: pixel -> degree
propixx_res = [1920 1080];                                   % propixx resolution
scr.w            = 72;                                       % screen width in cm
scr.h            = 40.5;                                     % screen height in cm
scr.d            = 142;
scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
scr.scrdegrw     = asind(scr.w/scr.ch);
scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pi

% load trigger
load(fullfile(pth,'experiment','trigdef.mat'))

% set up such that centre coordinates are centre of histogram bin
hist_bins_x = [sort(propixx_res(1)/2-scr.onedegrpix/4:-scr.onedegrpix/4:0),propixx_res(1)/2+scr.onedegrpix/4:scr.onedegrpix/4:propixx_res(1)];
hist_bins_y = [sort(propixx_res(2)/2-scr.onedegrpix/4:-scr.onedegrpix/4:0),propixx_res(2)/2+scr.onedegrpix/4:scr.onedegrpix/4:propixx_res(2)];

% saccade threshold
sac_thresh = 2;

s = 1;

% load trial structure
load(fullfile(trlpth,subj{s},"trl_overlap_meg_el_rsp.mat"))

% convert EDF to MAT and save as structure
file_in = fullfile(dtpth,subj{s},[subj{s}(end-3:end),'.edf']);
file_out = fullfile(dtpth,subj{s},'el_struct.mat');

el = kd_edf2mat2struct(edfconvpath,fullfile(dtpth,subj{s}),file_in,file_out);

% reject trials that could not be aligned with MEG
eltrl = elinfo.eltrl;                   % trial def
eltrl(elinfo.rej_rsp,:) = [];           % reject rsp specific trials 
elxy = elinfo.move;                     % xy over time
elxy(elinfo.rej_rsp,:,:) = [];            

eltrl = eltrl(~elinfo.rejtrl_all,:);
elxy = elxy(~elinfo.rejtrl_all,:,:);

% identify saccades
sacrej = [];
c = 0; % counter variable
for sac = 1:length(el.Events.Esacc.start)

    % starts within
    y = find((eltrl(:,3) < el.Events.Esacc.start(sac)) + (el.Events.Esacc.start(sac) < eltrl(:,4) - (0.4 *  elinfo.fsample)) == 2 );

    if ~isempty(y)
        c = c+1;

        % start,duration of saccade
        sac_start = el.Events.Ssacc.time(sac) - eltrl(y,2);
        sac_dur = el.Events.Esacc.duration(sac);

        sacrej(c,:) = y;
    end


end
sacrej = unique(sacrej);

% 2. find saccades >threshold
c = 0; % counter variable
sacrej_liberal = [];
for sac = 1:length(sacrej)

    % cut out eye position from end of baseline to end of trial
    bsl = 1.5* elinfo.fsample;                                % baseline dur
    trl_l = eltrl(sacrej(sac),end)-eltrl(sacrej(sac),end-2); % trial length
    cur_trl = squeeze(elxy(sacrej(sac),bsl:trl_l,:));
    % deviation x direction
    if ~isempty(find(cur_trl(:,1) < propixx_res(1)/2 - scr.onedegrpix*sac_thresh)) || ~isempty(find(cur_trl(:,1) > propixx_res(1)/2 + scr.onedegrpix*sac_thresh))
        c = c+1;

        sacrej_liberal(c,:) = sacrej(sac);

        % deviation y direction
    elseif ~isempty(find(cur_trl(:,2) < propixx_res(2)/2 - scr.onedegrpix*sac_thresh)) || ~isempty(find(cur_trl(:,2) > propixx_res(2)/2 + scr.onedegrpix*sac_thresh))
        c = c+1;

        sacrej_liberal(c,:) = sacrej(sac);
    end
end
sacrej_liberal = unique(sacrej_liberal);

%% Heatmap
clean_trl = setdiff(1:size(eltrl,1),sacrej_liberal);

pdf_clean = zeros(size(clean_trl,2),length(hist_bins_x)-1,length(hist_bins_y)-1);

% create probability density for each trial
for cltrl = 1:length(clean_trl)
    % extract eye position
    x_coord = elxy(clean_trl(cltrl),1.5* elinfo.fsample:end-0.4* elinfo.fsample,1);
    y_coord = elxy(clean_trl(cltrl),1.5* elinfo.fsample:end-0.4* elinfo.fsample,2);
    % discard zeros 
    x_coord = x_coord(x_coord ~= 0);
    y_coord = y_coord(x_coord ~= 0);
    % estimate probability denisty for this trial
    [pdf_clean(cltrl,:,:),Xedges,Yedges] = histcounts2(x_coord,y_coord,hist_bins_x,hist_bins_y,'Normalization','probability');
end

% repeat for saccade trials
pdf_sac = zeros(length(unique(sacrej_liberal)),length(hist_bins_x)-1,length(hist_bins_y)-1);

for sac = 1:length(sacrej_liberal)
    x_coord = elxy(sacrej_liberal(sac),1.5* elinfo.fsample:end-0.4* elinfo.fsample,1);
    y_coord = elxy(sacrej_liberal(sac),1.5* elinfo.fsample:end-0.4* elinfo.fsample,2);
    x_coord = x_coord(x_coord ~= 0);
    y_coord = y_coord(x_coord ~= 0);

    [pdf_sac(sac,:,:),Xedges,Yedges] = histcounts2(x_coord,y_coord,hist_bins_x,hist_bins_y,'Normalization','probability');
end


% colormap
cm = cbrewer('div','RdBu',21);
cm = cm(1:11,:);

m_pdf_clean = squeeze(mean(pdf_clean,1));
fig = figure;
subplot(211)
imagesc(Xedges./scr.onedegrpix, Yedges./scr.onedegrpix,m_pdf_clean)
hold on
% degree around fix dot
rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix)./scr.onedegrpix, 1*2, 1*2]);
hold on
% search display
rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix*5)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix*5)./scr.onedegrpix, 5*2, 5*2])

xticks(0:4:scr.scrdegrw)
yticks(0:3:scr.scrdegrh)
xlabel('x degree')
ylabel('y degree')
title(['no saccades ',num2str(length(clean_trl))])
colormap(flipud(cm))
caxis([0, round(max(m_pdf_clean(:)),1)])
cb = colorbar;
cb.Label.String = 'probability of coordinates';
axis xy
pbaspect([ceil(scr.scrdegrw)/ceil(scr.scrdegrh), 1,1])
if ~isempty(pdf_sac)
    subplot(212)
    m_pdf_sac = squeeze(mean(pdf_sac,1));
    imagesc(Xedges./scr.onedegrpix, Yedges./scr.onedegrpix,m_pdf_sac)
    % degree around fix dot
    rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix)./scr.onedegrpix, 1*2, 1*2]);
    hold on
    % search display
    rectangle('Position',[(propixx_res(1)/2-scr.onedegrpix*5)./scr.onedegrpix, (propixx_res(2)/2-scr.onedegrpix*5)./scr.onedegrpix, 5*2, 5*2])

    xlabel('x degree')
    ylabel('y degree')
    title(['saccades ',num2str(length(sacrej_liberal)),' trials'])
    xticks(0:4:scr.scrdegrw)
    yticks(0:3:scr.scrdegrh)
    colormap(flipud(cm))
    caxis([0, round(max(m_pdf_sac(:)),1)])
    cb = colorbar;
    cb.Label.String = 'probability of coordinates';
    axis xy
    pbaspect([ceil(scr.scrdegrw)/ceil(scr.scrdegrh), 1,1])
end
print(fig,fullfile(pltpth,[subj{s},'_heatmap_sac_thresh_',num2str(sac_thresh),'_degr']),'-dpng')
close all

% save to elinfo
elinfo.sactrial_lib = sacrej_liberal;
elinfo.pdf_clean = pdf_clean;
elinfo.pdf_sac = pdf_sac;

save(fullfile(trlpth,subj{s},"trl_overlap_meg_el_rsp.mat"),'elinfo','-append')
