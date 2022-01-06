%% VS + RFT
% PhD project 2

% inspect eye movement data

% [c] Katharina Duecker

clear all; close all; clc; beep off
% define paths
pth = 'E:\UoB\Proj2_Visual Search';
mtpth = 'Y:\Visual Search RFT\matlab scripts';
addpath(fullfile(mtpth,'edf-converter-master')); %edf2mat converter
% addpath('Y:\Visual Search RFT\matlab scripts\edfconverter')
addpath('C:\Program Files (x86)\SR Research\EyeLink\bin')
command = 'edf2asc.exe';
dtpth = fullfile(pth,'data');
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')            % fieldtrip
addpath('Y:\Visual Search RFT\matlab scripts\cbrewer')
pltpth = fullfile('Y:\Visual Search RFT\results\eyelink','figures');
elpth = fullfile('Y:\Visual Search RFT\results\eyelink','1 identify saccades');
prepropth = 'Y:\Visual Search RFT\results\meg\2 merged edf mat';

mkdir(pltpth)
mkdir(elpth)
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds
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
sac_thresh = [1.5,2,2.5];

num_trl_subj = zeros(length(subjfolds),length(sac_thresh));

for s = 1:length(subjfolds)
        trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
        elpth = fullfile('Y:\Visual Search RFT\results\eyelink','1 identify saccades');
        sac_thresh = [1,1.5,2,2.5,3,3.5];
        
        
        if exist(fullfile(trl_merge_pth, subjfolds{s},'trl_overlap_meg_el_rsp_sac_search.mat'))
            load(fullfile(trl_merge_pth, subjfolds{s},'trl_overlap_meg_el_rsp_sac_search.mat'))
        end
        
              
        for st = 1:length(sac_thresh)
            
            % 1. saccades
            sacrej = [];
            c = 0; % counter variable
            for sac = 1:length(el.Events.Esacc.end)
                
                % starts within
                y = find((eltrl(:,3) < el.Events.Esacc.start(sac)) + (el.Events.Esacc.start(sac) < eltrl(:,4) - (0.4 *  elinfo.fsample)) == 2 );
                
                if ~isempty(y)
                    c = c+1;
                    
                    % start,duration of saccade
                    sac_start = el.Events.Ssacc.time(sac) - eltrl(y,2);
                    sac_dur = el.Events.Esacc.duration(sac);
                    
                    sacrej(c,:) = [y,sac_start, sac_dur];
                end
                
                
            end
            
            
            % 2. find saccades >1 °
            c = 0; % counter variable
            sacrej_liberal = [];
            for sac = 1:length(sacrej)
                
                % loop over saccades - where does eye movement deviate from centre
                % by more than 1°?
                
                % cut out eye position from end of baseline to end of trial
                bsl = 1.5* elinfo.fsample;                                % baseline dur
                trl_l = eltrl(sacrej(sac),end)-eltrl(sacrej(sac),end-2); % trial length
                cur_trl = squeeze(elxy(sacrej(sac),bsl:trl_l,:));
                % deviation x direction
                if ~isempty(find(cur_trl(:,1) < propixx_res(1)/2 - scr.onedegrpix*sac_thresh(st))) || ~isempty(find(cur_trl(:,1) > propixx_res(1)/2 + scr.onedegrpix*sac_thresh(st)))
                    c = c+1;
                    
                    sacrej_liberal(c,:) = sacrej(sac);
                    
                    % deviation y direction
                elseif ~isempty(find(cur_trl(:,2) < propixx_res(2)/2 - scr.onedegrpix*sac_thresh(st))) || ~isempty(find(cur_trl(:,2) > propixx_res(2)/2 + scr.onedegrpix*sac_thresh(st)))
                    c = c+1;
                    
                    sacrej_liberal(c,:) = sacrej(sac);
                end
            end
            
            % blinks
            blinkrej = [];
            for sac = 1:length(el.Events.Sblink.time)
                x = find((eltrl(:,3) <= el.Events.Sblink.time(sac)) + (eltrl(:,4)-0.3* elinfo.fsample >= el.Events.Sblink.time(sac)) == 2 );
                
                blinkrej = [blinkrej;x];
            end
            blinkrej = unique(blinkrej);
            
            %% How much time do people spend at fix dot?
            
            % probability density clean trials (baseline vs during search)
            clean_trl = setdiff(1:size(eltrl,1),unique([sacrej_liberal;blinkrej]));
            
            pdf_clean = zeros(size(clean_trl,2),length(hist_bins_x)-1,length(hist_bins_y)-1);
            
            % create probability density for each trial
            for cltrl = 1:length(clean_trl)
                % extract eye position
                x_coord = elxy(clean_trl(cltrl),1.5* elinfo.fsample:end-0.4* elinfo.fsample,1);
                y_coord = elxy(clean_trl(cltrl),1.5* elinfo.fsample:end-0.4* elinfo.fsample,2);
                % discard zeros (duration up to button press)
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
            
            print(fig,fullfile(pltpth,[subjfolds{s},'_heatmap_Esacc_sac_thresh_',num2str(sac_thresh(st)*10),'_degr']),'-dpng')
            
         
            elinfo.sactrial_lib = sacrej_liberal;
            elinfo.blink = blinkrej;
            elinfo.pdf_clean = pdf_clean;
            elinfo.pdf_sac = pdf_sac;
            num_trl_subj(s,st) = length(clean_trl);
            
            
            save(fullfile(elpth,[subjfolds{s},'_elinfo_sac_thresh_',num2str(sac_thresh(st)*10),'.mat']),'elinfo')
        end
        clearvars  N* b* a* x* y* c* e* f* X* Y* idx* d trl* sac*
        close all
end

elpth = fullfile('Y:\Visual Search RFT\results\eyelink','1 identify saccades');

save(fullfile(elpth,'num_trl_thresh.mat'),'num_trl_subj')