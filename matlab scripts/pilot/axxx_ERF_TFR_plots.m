%% Pilot analysis: VS, alpha and RFT
% plots TFR

%% Settings
clear all; close all; clc; beep off;
mpth = 'X:\';
cd(mpth)
megpth = fullfile(mpth,'pilot','meg');
rsppth = fullfile(mpth,'pilot','responses');
cdpth = fullfile(mpth,'experiment');
mtlpth = fullfile(mpth,'matlab','pilot');
addpath(mtlpth)
addpath('C:\Users\katha\Documents\MATLAB\fieldtrip')
ft_defaults;
subjcode = 'b57a';
pthout = fullfile(mpth,'pilot','results',subjcode);
pthfig = fullfile(mpth,'pilot','results','plots',subjcode);
load(fullfile(pthout,'TFR_ERF.mat'))

% plot SOI6067

cfg = [];
cfg.baseline = [-1.25 -.25];
cfg.baselinetype = 'relchange';
cfg.channel = SOI6067;
cfg.ylim = [55 75];
cfg.colorbar = 'yes';
ft_singleplotTFR(cfg,TFR)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
cb = colorbar;
cb.Label.String = 'power relative to baseline (a.u.)';
hold on
line(TFR.time,repmat(60,1,length(TFR.time)),'Color','black')
hold on
line(TFR.time,repmat(67,1,length(TFR.time)),'Color','black')
xlim([-1 0.5])

print(fullfile(pthfig,'SOI6067'),'-dpng')
close all


for s = 1:length(SOI60)
    subplot(3,2,s)
    cfg.channel = SOI60{s};
    ft_singleplotTFR(cfg,TFR)
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
    cb = colorbar;
    cb.Label.String = 'power relative to baseline (a.u.)';
    hold on
    line(TFR.time,repmat(60,1,length(TFR.time)),'Color','black','LineWidth',1)
    hold on
    line(TFR.time,repmat(67,1,length(TFR.time)),'Color','black','LineWidth',1)
    xlim([-1 0.5])
end

print(fullfile(pthfig,'SOI60'),'-dpng')
close all

h = 0;
fig = figure;
set(fig,'Position',[0 0 1550 800])
for s = 1:length(SOI67)
    x = mod(s,6);
    if x == 0
        x = 6;
    end
    subplot(3,2,x)
    cfg.channel = SOI67{s};
    ft_singleplotTFR(cfg,TFR)
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
    cb = colorbar;
    cb.Label.String = 'power relative to baseline (a.u.)';
    hold on
    line(TFR.time,repmat(60,1,length(TFR.time)),'Color','black','LineWidth',1)
    hold on
    line(TFR.time,repmat(67,1,length(TFR.time)),'Color','black','LineWidth',1)
    xlim([-1 0.5])
    if ~mod(s,6) || s == length(SOI67)
        h = h+1;
        print(fullfile(pthfig,['SOI67_',num2str(h)]),'-dpng')
            pause

        close all
        fig = figure;
        set(fig,'Position',[0 0 1550 800])
    end
end


%% Topoplots

fig = figure;
set(fig,'Position',[0 0 1550 800])
% 60
cfg = [];
cfg.xlim = [0.1 0.5];
cfg.ylim = [59 60];
cfg.baseline = [-1.25 -.25];
cfg.baselinetype = 'relchange';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.colorbar = 'yes';

subplot(121)

ft_topoplotTFR(cfg,TFR)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
cb = colorbar;
cb.Label.String = 'power relative to baseline (a.u.)';

% 67
subplot(122)
cfg.ylim = [66 67];
ft_topoplotTFR(cfg,TFR)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
cb = colorbar;
cb.Label.String = 'power relative to baseline (a.u.)';

print(fig,fullfile(pthfig,'topo_RFT'),'-dpng')
close all