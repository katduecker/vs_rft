%% VS + RFT
% PhD project 2

% Coherence contrasts
clear all; close all; clc; beep off;
rej_sac = 0;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
addpath(fullfile(pth,'matlab scripts/','cbrewer'))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh');

cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi');
mkdir(cohfigpth)
maxfpth = fullfile(pth,'results','meg', '1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;


s = 1;
condi = {{'6067','ni'},{'6067','ti'}};                 % conditions to be contrasted

toi = [-2.5 2];                             % start and end of trial in sec

fw = 2;                                     % bandwidth bp filter
fs = 1000;

if length(condi{1}) > 1
    Tfreq = [str2double(condi{1}{1}(1:2)) str2double(condi{2}{1}(1:2))];
    Dfreq = [str2double(condi{1}{1}(3:4)) str2double(condi{2}{1}(3:4))];
else
    
    Tfreq = [str2double(condi{1}{1}(1:2)) str2double(condi{2}{1}(1:2))];
    Dfreq = [str2double(condi{1}{1}(3:4)) str2double(condi{2}{1}(3:4))];
    
end



%plot and contrast coherence
for r = 1:length(rej_sac)
    for fwc = 1:length(fw)
        
        if rej_sac(r)
            cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','rej_sac');
        else
            cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh');
        end
            
        d = dir(cohpth);
        d = {d.name};
        subj = d(strncmp(d,'202',3));
        %% create place holder struct
        

        
        
        % list maxfiltered data
        d = dir(fullfile(maxfpth,subj{s}));
        f = {d.name};
        % find fif file
        idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
        idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
        f = f(idxx);
        
        

        
        % trial structure to load in trl
        load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))
        
        trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
        trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;
        
        cfg = [];
        cfg.dataset = fullfile(maxfpth,subj{s},f{1});
        cfg.preproc.detrend = 'yes';
        cfg.trl = trlstruct{1};
        cfg.channel = {'MEG','MISC004','MISC005'};
        % load in data for this part
        data = ft_preprocessing(cfg);
        
        %% Placeholder coherence
        
        % fourier transform
        cfg            = [];
        cfg.output     = 'pow';
        cfg.method     = 'mtmconvol';
        %cfg.pad        = 'nextpow2';
        cfg.taper      = 'hanning';
        cfg.toi        = -1.5:2;
        cfg.foi        = 53:78;
        cfg.t_ftimwin  = ones(length(cfg.foi),1).*.1;
        cfg.tapsmofrq  = 3;
        cfg.keeptrials = 'no';
        cfg.channel    = {'MEGGRAD', 'MISC004','MISC005'};
        freqgrad       = ft_freqanalysis(cfg, data);
        cfg.channel    = {'MEGMAG', 'MISC004'};
        freqmag        = ft_freqanalysis(cfg, data);
        freqgrad.time  = toi(1):1/fs:toi(2);
        freqgrad.freq  = 53:78;
        freqmag.time  = toi(1):1/fs:toi(2);
        freqmag.freq  = 53:78;
        
        if rej_sac(r)
            respth = fullfile(cohfigpth,'rej_sac',subj{s},['fw ',num2str(fw(fwc)),' Hz'], [strjoin(condi{1}),' vs ',strjoin(condi{2})]);
        else
            respth = fullfile(cohfigpth,subj{s},['fw ',num2str(fw(fwc)),' Hz'], [strjoin(condi{1}),' vs ',strjoin(condi{2})]);
        end
        
        mkdir(respth)
        
        %% Load Coherence
        
        % load conditions
        for c = 1:length(condi)
            curcond = ['coh_',strjoin(condi{c},'_'),'_freqw_',num2str(fw(fwc))];
            dat{c} = load(fullfile(cohpth,subj{s},curcond));
        end
        
        % load soi
        
        load(fullfile(soipth,subj{s},'soi_stat.mat'))
        
        soimag = soi_stat(logical(cell2mat(cellfun(@(x) strcmp(x(end),'1'),soi_stat,'UniformOutput',false))));
        soigrad = soi_stat(~ismember(soi_stat,soimag));
        
        
        % average display onset to minimum RT     
        minrt = min(cell2mat(rspinfo.trl(:,3)));
        
        
        
        %% GRADIOMETERS
        % plot
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.layout = 'neuromag306planar.lay';
        cfg.channel = 'MEGGRAD';
        cfg.marker = 'off';
        cfg.highlight = 'labels';
        cfg.highlightchannel = soigrad;
        cfg.highlightsymbol = 'o';
        cfg.highlightsize = 4;
        cfg.highlightcolor = [0 0.5 0];
        cfg.zlim = 'zeromax';
        cfg.xlim = [0.1 minrt];
        cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
        
        
        % T Frequency
        freqgrad.powspctrm = dat{1}.coh.cohTgrad;
        cfg.ylim = [Tfreq(1) Tfreq(1)];
        cfg.colorbar = 'yes';
                cfg.colorbartext = 'coh';

        ft_topoplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Tdiode_grad_',strjoin(condi{1}),'.png']))
        % cfg.ylim = [67 67];
        % ft_topoplotTFR(cfg,freqgrad)
        % set(gcf,'Position',[0 0 1920 1080])
        %
        % saveas(gcf,fullfile(respth,'T60Hz_coh_Tdiode_67Hz_grad_sanity.png'))
        close all
        
        % D Frequency
        freqgrad.powspctrm = dat{1}.coh.cohDgrad;
        cfg.ylim = [Dfreq(1) Dfreq(1)];
        cfg.zlim = 'zeromax';
        ft_topoplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Ddiode_grad_',strjoin(condi{1}),'.png']))
        
        % cfg.ylim = [60 60];
        % cfg.zlim = 'zeromax';
        % ft_topoplotTFR(cfg,freqgrad)
        % set(gcf,'Position',[0 0 1920 1080])
        %
        % saveas(gcf,fullfile(respth,'D67Hz_coh_Ddiode_60Hz_grad_sanity.png'))
        
        close all
        
        % condition 2
        % T Frequency
        freqgrad.powspctrm = dat{2}.coh.cohTgrad;
        % plot
        cfg.ylim = [Tfreq(2) Tfreq(2)];
        ft_topoplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Tdiode_grad_',strjoin(condi{2}),'.png']))
        
        % cfg.ylim = [60 60];
        % ft_topoplotTFR(cfg,freqgrad)
        % set(gcf,'Position',[0 0 1920 1080])
        %
        % saveas(gcf,fullfile(respth,'T67Hz_coh_Tdiode_60Hz_grad_sanity.png'))
        %
        % close all
        
        % D 60 Hz
        freqgrad.powspctrm = dat{2}.coh.cohDgrad;
        cfg.ylim = [Dfreq(2) Dfreq(2)];
        ft_topoplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Ddiode_grad_',strjoin(condi{2}),'.png']))
        
        % % saveas(gcf,fullfile(respth,'D60Hz_coh_Ddiode_grad.png'))
        % cfg.ylim = [67 67];
        % ft_topoplotTFR(cfg,freqgrad)
        % set(gcf,'Position',[0 0 1920 1080])
        %
        % saveas(gcf,fullfile(respth,'D60Hz_coh_Ddiode_67Hz_grad_sanity.png'))
        
        %% Singleplot
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.layout = 'neuromag306planar.lay';
        cfg.channel = soigrad;
        cfg.title = 'T unknown 60 Hz';
        cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
        cfg.colorbar = 'yes';
        cfg.colorbartext = 'coh';
        freqgrad.powspctrm = dat{1}.coh.cohTgrad;
        ft_singleplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['TFRcoh_T_grad_',strjoin(condi{1},'_'),'.png']))
        close all
        cfg.title = 'T known 60 Hz';
        freqgrad.powspctrm = dat{2}.coh.cohTgrad;
        ft_singleplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['TFRcoh_T_grad_',strjoin(condi{2},'_'),'.png']))
        close all
        cfg.title = 'T unknown vs known 60 Hz';
        cfg.colorbartext = '(coh unknwon - coh nown)/(coh unknown + coh known)';
        cfg.colormap = cm;
        freqgrad.powspctrm = (round(dat{1}.coh.cohTgrad,3)-round(dat{2}.coh.cohTgrad,3))./(round(dat{1}.coh.cohTgrad,3)+round(dat{2}.coh.cohTgrad,3));
        ft_singleplotTFR(cfg,freqgrad)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['TFRcontrast_T_grad_',strjoin(condi{1},'_'),'_vs_',strjoin(condi{2},'_'),'.png']))
        
        %% Contrast conditions
        

        avgsamp = abs(toi(1))*fs:abs(toi(1))*fs+round(minrt,3)*fs;
        
        for d = 1:length(dat)
            
            avg{d}.cohTgrad = mean(dat{d}.coh.cohTgrad(:,:,avgsamp),3);
            avg{d}.cohDgrad = mean(dat{d}.coh.cohDgrad(:,:,avgsamp),3);
            avg{d}.cohTmag = mean(dat{d}.coh.cohTmag(:,:,avgsamp),3);
            avg{d}.cohDmag = mean(dat{d}.coh.cohDmag(:,:,avgsamp),3);
        end
        
        
        if Tfreq(1) == Dfreq(2)
            
            freqgrad.powspctrm = (round(avg{1}.cohTgrad,3)-round(avg{2}.cohDgrad,3))./(round(avg{1}.cohTgrad,3)+round(avg{2}.cohDgrad,3));
            freqgrad.time = 0;
            
            cfg.xlim = [];
            cfg.ylim = [Tfreq(1) Tfreq(1)];
            cfg.zlim = [-max(abs(freqgrad.powspctrm(:))) max(abs(freqgrad.powspctrm(:)))];
            cfg.colorbar = 'yes';
            cfg.colorbartext = ['coh T vs D ',num2str(Tfreq(1)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqgrad)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_TvsD_',strjoin(condi{1},'_'),'_vs_',strjoin(condi{2},'_'),'_grad.png']))
            close all
            freqgrad.powspctrm = (round(avg{2}.cohTgrad,3)-round(avg{1}.cohDgrad,3))./(round(avg{2}.cohTgrad,3)+round(avg{1}.cohDgrad,3));
            freqgrad.time = 0;
            
            cfg.xlim = [];
            cfg.ylim = [Tfreq(2) Tfreq(2)];
            cfg.zlim = [-max(abs(freqgrad.powspctrm(:))) max(abs(freqgrad.powspctrm(:)))];
            cfg.colorbartext = ['coh T vs D ',num2str(Tfreq(2)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqgrad)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_TvsD_',strjoin(condi{2},'_'),'_vs_',strjoin(condi{1},'_'),'_grad.png']))
            close all
            
        elseif Tfreq(1) == Tfreq(2)
            
            freqgrad.powspctrm = (round(avg{1}.cohTgrad,3)-round(avg{2}.cohTgrad,3))./(round(avg{1}.cohTgrad,3)+round(avg{2}.cohTgrad,3));
            freqgrad.time = 0;
            
            cfg.xlim = [];
            cfg.ylim = [Tfreq(1) Tfreq(1)];
            cfg.zlim = [-max(abs(freqgrad.powspctrm(:))) max(abs(freqgrad.powspctrm(:)))];
            cfg.colorbar = 'yes';
            cfg.colorbartext = ['coh T known vs unknown',num2str(Tfreq(1)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqgrad)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_Tknownvsunknown_',strjoin(condi{1},'_'),'_vs_',strjoin(condi{2},'_'),'_grad.png']))
            close all
            freqgrad.powspctrm = (round(avg{1}.cohDgrad,3)-round(avg{2}.cohDgrad,3))./(round(avg{1}.cohDgrad,3)+round(avg{2}.cohDgrad,3));
            freqgrad.time = 0;
            
            cfg.xlim = [];
            cfg.ylim = [Tfreq(2) Tfreq(2)];
            cfg.zlim = [-max(abs(freqgrad.powspctrm(:))) max(abs(freqgrad.powspctrm(:)))];
            cfg.colorbartext = ['coh D known vs unknown ',num2str(Dfreq(1)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqgrad)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_Dknownvsunknown_',strjoin(condi{2},'_'),'_vs_',strjoin(condi{1},'_'),'_grad.png']))
            close all
        end
        
        
        % grad
        figure('Position',[0,0,1920,1080])
        for g = 1:length(soigrad)-1
            subplot(4,3,g)
            l = strcmp(freqgrad.label,soigrad{g});
            
            plot(freqgrad.freq,avg{1}.cohTgrad(l,:),'Color',[0 0.4470 0.7410],'LineWidth',1)
            hold on
            plot(freqgrad.freq,avg{1}.cohDgrad(l,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
            hold on
            plot(freqgrad.freq,avg{2}.cohTgrad(l,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
            hold on
            plot(freqgrad.freq,avg{2}.cohDgrad(l,:),'Color',[0.6350 0.0780 0.1840],'LineWidth',1)
            
            plot(Tfreq(1),avg{1}.cohTgrad(l,freqgrad.freq==Tfreq(1)),'*','Color',[0 0.4470 0.7410])
            plot(Dfreq(1),avg{1}.cohDgrad(l,freqgrad.freq==Dfreq(1)),'*','Color',[0.8500 0.3250 0.0980])
            plot(Tfreq(2),avg{2}.cohTgrad(l,freqgrad.freq==Tfreq(2)),'*','Color',[0.4660 0.6740 0.1880])
            plot(Dfreq(2),avg{2}.cohDgrad(l,freqgrad.freq==Dfreq(2)),'*','Color',[0.6350 0.0780 0.1840])
            title(soigrad{g})
        end
        g = g+1;
        subplot(4,3,g:g+1)
        l = strcmp(freqgrad.label,soigrad{g});
        plot(freqgrad.freq,avg{1}.cohTgrad(l,:),'Color',[0 0.4470 0.7410],'LineWidth',1)
        hold on
        plot(freqgrad.freq,avg{1}.cohDgrad(l,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold on
        plot(freqgrad.freq,avg{2}.cohTgrad(l,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
        hold on
        plot(freqgrad.freq,avg{2}.cohDgrad(l,:),'Color',[0.6350 0.0780 0.1840],'LineWidth',1)
        
        plot(Tfreq(1),avg{1}.cohTgrad(l,freqgrad.freq==Tfreq(1)),'*','Color',[0 0.4470 0.7410])
        plot(Dfreq(1),avg{1}.cohDgrad(l,freqgrad.freq==Dfreq(1)),'*','Color',[0.8500 0.3250 0.0980])
        plot(Tfreq(2),avg{2}.cohTgrad(l,freqgrad.freq==Tfreq(2)),'*','Color',[0.4660 0.6740 0.1880])
        plot(Dfreq(2),avg{2}.cohDgrad(l,freqgrad.freq==Dfreq(2)),'*','Color',[0.6350 0.0780 0.1840])
        title(soigrad{g})
        if length(condi{1}) > 1
            legend({['T',num2str(Tfreq(1)),' Hz known'],['D',num2str(Dfreq(1)),' Hz known'],['T',num2str(Tfreq(2)),' Hz unknown'],['D',num2str(Dfreq),' Hz unknown']},'AutoUpdate','off','Location','bestoutside')
        else
            legend({['T',num2str(Tfreq(1)),' Hz'],['D',num2str(Dfreq(1)),' Hz'],['T',num2str(Tfreq(2)),' Hz'],['D',num2str(Dfreq(2)),' Hz']},'AutoUpdate','off','Location','bestoutside')
            
        end
        saveas(gcf,fullfile(respth,['avg_coh_grad_',strjoin(condi{1}),'_vs_',strjoin(condi{2}),'.png']))
        
        close all
        
        %% MAGNETOMETERS
        % plot
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.layout = 'neuromag306mag.lay';
        cfg.channel = 'MEGGRAD';
        cfg.marker = 'off';
        cfg.highlight = 'labels';
        cfg.highlightchannel = soimag;
        cfg.highlightsymbol = 'o';
        cfg.highlightsize = 4;
        cfg.highlightcolor = [0 0.5 0];
        cfg.zlim = 'zeromax';
        cfg.xlim = [0.1 minrt];
        cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
        cfg.colorbar = 'yes';
        
        % T Frequency
        freqmag.powspctrm = dat{1}.coh.cohTmag;
        cfg.ylim = [Tfreq(1) Tfreq(1)];
        cfg.colorbar = 'yes';
        cfg.colorbartext = 'coh';
        ft_topoplotTFR(cfg,freqmag)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Tdiode_mag_',strjoin(condi{1}),'.png']))
        % cfg.ylim = [67 67];
        % ft_topoplotTFR(cfg,freqgrad)
        % set(gcf,'Position',[0 0 1920 1080])
        %
        % saveas(gcf,fullfile(respth,'T60Hz_coh_Tdiode_67Hz_grad_sanity.png'))
        close all
        
        % D Frequency
        freqmag.powspctrm = dat{1}.coh.cohDmag;
        cfg.ylim = [Dfreq(1) Dfreq(1)];
        cfg.zlim = 'zeromax';
        ft_topoplotTFR(cfg,freqmag)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Ddiode_mag_',strjoin(condi{1}),'.png']))
        
        
        close all
        
        % condition 2
        % T Frequency
        freqmag.powspctrm = dat{2}.coh.cohTmag;
        % plot
        cfg.ylim = [Tfreq(2) Tfreq(2)];
        ft_topoplotTFR(cfg,freqmag)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Tdiode_mag_',strjoin(condi{2}),'.png']))
        
        close all
        
        % D 60 Hz
        freqmag.powspctrm = dat{2}.coh.cohDmag;
        cfg.ylim = [Dfreq(2) Dfreq(2)];
        ft_topoplotTFR(cfg,freqmag)
        set(gcf,'Position',[0 0 1920 1080])
        saveas(gcf,fullfile(respth,['Ddiode_mag_',strjoin(condi{2}),'.png']))
        
        close all
        %% Contrast conditions
        
        if Tfreq(1) == Tfreq(2)
            
            freqmag.powspctrm = (round(avg{1}.cohTmag,3)-round(avg{2}.cohTmag,3))./(round(avg{1}.cohTmag,3)+round(avg{2}.cohTmag,3));
            %freqmag.powspctrm = (round(avg{1}.cohTmag,3)./round(avg{2}.cohTmag,3)) - 1;
            freqmag.time = 0;
            
            cfg.xlim = [];
            cfg.ylim = [Tfreq(1) Tfreq(1)];
            cfg.zlim = [-max(abs(freqmag.powspctrm(:))) max(abs(freqmag.powspctrm(:)))];
            cfg.colorbartext = ['coh T known vs unknown ',num2str(Tfreq(1)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqmag)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_T known vs unknown_',strjoin(condi{1},'_'),'_vs_',strjoin(condi{2},'_'),'_mag.png']))
            close all
            
            freqmag.powspctrm = (round(avg{1}.cohDmag,3)-round(avg{2}.cohDmag,3))./(round(avg{1}.cohDmag,3)+round(avg{2}.cohDmag,3));
            % freqmag.powspctrm = (round(avg{2}.cohTmag,3)./round(avg{1}.cohTmag,3)-1;
            cfg.xlim = [];
            cfg.ylim = [Tfreq(2) Tfreq(2)];
            cfg.zlim = [-max(abs(freqmag.powspctrm(:))) max(abs(freqmag.powspctrm(:)))];
            
            cfg.colorbartext = ['coh D known vs unknown',num2str(Dfreq(1)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqmag)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_D known vs unknown',strjoin(condi{2},'_'),'_vs_',strjoin(condi{1},'_'),'_mag.png']))
            close all
            
        elseif Tfreq(1) == Dfreq(2)
            
            freqmag.powspctrm = (round(avg{1}.cohTmag,3)-round(avg{2}.cohDmag,3))./(round(avg{1}.cohTmag,3)+round(avg{2}.cohDmag,3));
            %freqmag.powspctrm = (round(avg{1}.cohTmag,3)./round(avg{2}.cohTmag,3)) - 1;
            freqmag.time = 0;
            
            cfg.xlim = [];
            cfg.ylim = [Tfreq(1) Tfreq(1)];
            cfg.zlim = [-max(abs(freqmag.powspctrm(:))) max(abs(freqmag.powspctrm(:)))];
            cfg.colorbartext = ['coh T vs D ',num2str(Tfreq(1)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqmag)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_TvsD_',strjoin(condi{1},'_'),'_vs_',strjoin(condi{2},'_'),'_mag.png']))
            close all
            
            freqmag.powspctrm = (round(avg{1}.cohDmag,3)-round(avg{2}.cohDmag,3))./(round(avg{1}.cohDmag,3)+round(avg{2}.cohDmag,3));
            % freqmag.powspctrm = (round(avg{2}.cohTmag,3)./round(avg{1}.cohTmag,3)-1;
            cfg.xlim = [];
            cfg.ylim = [Tfreq(2) Tfreq(2)];
            cfg.zlim = [-max(abs(freqmag.powspctrm(:))) max(abs(freqmag.powspctrm(:)))];
            
            cfg.colorbartext = ['coh T vs D',num2str(Tfreq(2)),' Hz'];
            cfg.colormap = cm;
            
            ft_topoplotTFR(cfg,freqmag)
            set(gcf,'Position',[0 0 1920 1080])
            
            saveas(gcf,fullfile(respth,['relcoh_DvsT',strjoin(condi{2},'_'),'_vs_',strjoin(condi{1},'_'),'_mag.png']))
            close all
        end
        
        
        figure('Position',[0,0,1920,1080])
        for g = 1:length(soimag)-1
            subplot(4,3,g)
            l = strcmp(freqmag.label,soimag{g});
            
            plot(freqmag.freq,avg{1}.cohTmag(l,:),'Color',[0 0.4470 0.7410],'LineWidth',1)
            hold on
            plot(freqmag.freq,avg{1}.cohDmag(l,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
            hold on
            plot(freqmag.freq,avg{2}.cohTmag(l,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
            hold on
            plot(freqmag.freq,avg{2}.cohDmag(l,:),'Color',[0.6350 0.0780 0.1840],'LineWidth',1)
            
            plot(Tfreq(1),avg{1}.cohTmag(l,freqmag.freq==Tfreq(1)),'*','Color',[0 0.4470 0.7410])
            plot(Dfreq(1),avg{1}.cohDmag(l,freqmag.freq==Dfreq(1)),'*','Color',[0.8500 0.3250 0.0980])
            plot(Tfreq(2),avg{2}.cohTmag(l,freqmag.freq==Tfreq(2)),'*','Color',[0.4660 0.6740 0.1880])
            plot(Dfreq(2),avg{2}.cohDmag(l,freqmag.freq==Dfreq(2)),'*','Color',[0.6350 0.0780 0.1840])
            title(soimag{g})
        end
        g = g+1;
        subplot(4,3,g:g+1)
        l = strcmp(freqmag.label,soimag{g});
        plot(freqmag.freq,avg{1}.cohTmag(l,:),'Color',[0 0.4470 0.7410],'LineWidth',1)
        hold on
        plot(freqmag.freq,avg{1}.cohDmag(l,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
        hold on
        plot(freqmag.freq,avg{2}.cohTmag(l,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
        hold on
        plot(freqmag.freq,avg{2}.cohDmag(l,:),'Color',[0.6350 0.0780 0.1840],'LineWidth',1)
        
        plot(Tfreq(1),avg{1}.cohTmag(l,freqmag.freq==Tfreq(1)),'*','Color',[0 0.4470 0.7410])
        plot(Dfreq(1),avg{1}.cohDmag(l,freqmag.freq==Dfreq(1)),'*','Color',[0.8500 0.3250 0.0980])
        plot(Tfreq(2),avg{2}.cohTmag(l,freqmag.freq==Tfreq(2)),'*','Color',[0.4660 0.6740 0.1880])
        plot(Dfreq(2),avg{2}.cohDmag(l,freqmag.freq==Dfreq(2)),'*','Color',[0.6350 0.0780 0.1840])
        title(soimag{g})
        if length(condi{1}) > 1
            legend({['T',num2str(Tfreq(1)),' Hz known'],['D',num2str(Dfreq(1)),' Hz known'],['T',num2str(Tfreq(2)),' Hz unknown'],['D',num2str(Dfreq),' Hz unknown']},'AutoUpdate','off','Location','bestoutside')
        else
            legend({['T',num2str(Tfreq(1)),' Hz'],['D',num2str(Dfreq(1)),' Hz'],['T',num2str(Tfreq(2)),' Hz'],['D',num2str(Dfreq(2)),' Hz']},'AutoUpdate','off','Location','bestoutside')
            
        end
        saveas(gcf,fullfile(respth,['avg_coh_mag_',strjoin(condi{1}),'_vs_',strjoin(condi{2}),'.png']))
        close all
    end
end
%

% 
% saveas(gcf,fullfile(respth,['relative_coh_change_T60Hz_known_vs_unknown.png']))
% 
% %% MAGNETOMETERS
% 
% % plot
% cfg = [];
% cfg.parameter = 'powspctrm';
% cfg.layout = 'neuromag306mag.lay';
% cfg.zlim = 'zeromax';
% cfg.xlim = [0.1 1.1];
% cfg.colorbartext = 'coh';
% cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
% cfg.marker = 'off';
% cfg.highlight = 'labels';
% cfg.highlightchannel = soimag;
% cfg.highlightsymbol = 'o';
% cfg.highlightsize = 4;
% cfg.highlightcolor = [0 0.5 0];
% 
% % T 60 H
% freqmag.powspctrm = dat{1}.coh.cohTmag;
% cfg.ylim = [60 60];
% cfg.colorbar = 'yes';
% ft_topoplotTFR(cfg,freqmag)
% set(gcf,'Position',[0 0 1920 1080])
% 
% saveas(gcf,fullfile(respth,'T60Hz_coh_Tdiode_mag.png'))
% 
% close all
% % D 67 Hz
% freqmag.powspctrm = dat{1}.coh.cohDmag;
% cfg.ylim = [67 67];
% ft_topoplotTFR(cfg,freqmag)
% set(gcf,'Position',[0 0 1920 1080])
% 
% saveas(gcf,fullfile(respth,'D67Hz_coh_Ddiode_mag.png'))
% 
% % T 67 H
% freqmag.powspctrm = dat{2}.coh.cohTmag;
% cfg.ylim = [67 67];
% cfg.colorbar = 'yes';
% ft_topoplotTFR(cfg,freqmag)
% set(gcf,'Position',[0 0 1920 1080])
% 
% saveas(gcf,fullfile(respth,'T67Hz_coh_Tdiode_mag.png'))
% 
% close all
% % D 67 Hz
% freqmag.powspctrm = dat{2}.coh.cohDmag;
% cfg.ylim = [60 60];
% ft_topoplotTFR(cfg,freqmag)
% set(gcf,'Position',[0 0 1920 1080])
% 
% saveas(gcf,fullfile(respth,'D60Hz_coh_Ddiode_mag.png'))
% 
% %% Explore coherence over time
% 
% % GRAD
% 
% % 60 Hz
% freqgrad.powspctrm = dat{1}.coh.cohTgrad;
% 
% figure('Position',[0,0,1920,1080])
% for g = 1:length(soigrad)
%     subplot(4,3,g)
%     l = strcmp(freqgrad.label,soigrad{g});
%     imagesc(freqgrad.time,freqgrad.freq,squeeze(freqgrad.powspctrm(l,:,:)))
%     axis xy
%     colormap(cm(floor(length(cm)/2)+1:end,:))
%     colorbar
%     %caxis([-0.5 0.5])
%     title(soigrad{g})
% end
% suptitle('Target 60 Hz')
% 
% saveas(gcf,fullfile(respth,'T60Hz_grad_coh_time.png'))
% 
% % 67 Hz
% freqgrad.powspctrm = dat{1}.coh.cohDgrad;
% 
% figure('Position',[0,0,1920,1080])
% for g = 1:length(soigrad)
%     subplot(4,3,g)
%     l = strcmp(freqgrad.label,soigrad{g});
%     imagesc(freqgrad.time,freqgrad.freq,squeeze(freqgrad.powspctrm(l,:,:)))
%     axis xy
%     colormap(cm(floor(length(cm)/2)+1:end,:))
%     colorbar
%     %caxis([-0.5 0.5])
%     title(soigrad{g})
% end
% suptitle('Distractor 67 Hz')
% 
% saveas(gcf,fullfile(respth,'D67Hz_grad_coh_time.png'))
% 
% % MAG
% % 60 Hz
% freqmag.powspctrm = dat{1}.coh.cohTmag;
% 
% figure('Position',[0,0,1920,1080])
% for g = 1:length(soimag)
%     subplot(2,3,g)
%     l = strcmp(freqmag.label,soimag{g});
%     imagesc(freqmag.time,freqmag.freq,squeeze(freqmag.powspctrm(l,:,:)))
%     axis xy
%     colormap(cm(floor(length(cm)/2)+1:end,:))
%     colorbar
%     %caxis([-0.5 0.5])
%     title(soimag{g})
% end
% suptitle('Target 60 Hz')
% 
% saveas(gcf,fullfile(respth,'T60Hz_mag_coh_time.png'))
% 
% % 67 Hz
% freqmag.powspctrm = dat{1}.coh.cohDmag;
% 
% figure('Position',[0,0,1920,1080])
% for g = 1:length(soimag)
%     subplot(2,3,g)
%     l = strcmp(freqmag.label,soimag{g});
%     imagesc(freqmag.time,freqmag.freq,squeeze(freqmag.powspctrm(l,:,:)))
%     axis xy
%     colormap(cm(floor(length(cm)/2)+1:end,:))
%     colorbar
%     %caxis([-0.5 0.5])
%     title(soimag{g})
% end
% suptitle('Distractor 67 Hz')
% 
% saveas(gcf,fullfile(respth,'D67Hz_mag_coh_time.png'))

%% Average 



% % mag
% figure('Position',[0,0,1920,540])
% for g = 1:length(soimag)
%     subplot(2,3,g)
%     l = strcmp(freqmag.label,soimag{g});
%     avgcohT1 = mean(squeeze(dat{1}.coh.cohTmag(l,:,1600:2600)),2);
%     avgcohD1 = mean(squeeze(dat{1}.coh.cohDmag(l,:,1600:2600)),2);
%     avgcohT2 = mean(squeeze(dat{2}.coh.cohTmag(l,:,1600:2600)),2);
%     avgcohD2 = mean(squeeze(dat{2}.coh.cohDmag(l,:,1600:2600)),2);
%     
%     plot(freqgrad.freq,avgcohT1,'Color',[0 0.4470 0.7410],'LineWidth',1)
%     hold on 
%     plot(freqgrad.freq,avgcohD1,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
%     hold on
%     plot(freqgrad.freq,avgcohT2,'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
%     hold on 
%     plot(freqgrad.freq,avgcohD2,'Color',[0.6350 0.0780 0.1840],'LineWidth',1)
%     legend({'T60 Hz','D67 Hz','T67 Hz','D60 Hz'},'AutoUpdate','off')
%     
%     plot(60,avgcohT1(freqgrad.freq==60),'*','Color',[0 0.4470 0.7410])
%     plot(67,avgcohD1(freqgrad.freq==67),'*','Color',[0.8500 0.3250 0.0980])
%     plot(67,avgcohT2(freqgrad.freq==67),'*','Color',[0.4660 0.6740 0.1880])
%     plot(60,avgcohD2(freqgrad.freq==60),'*','Color',[0.6350 0.0780 0.1840])
%     
%     title(soimag{g})
% end
% saveas(gcf,fullfile(respth,'avg_coh_mag.png'))



%% cluster based statistics
