%% Make plot for RMS envelope as function of CVE (and related plots)
% Compute Coefficient of Variation (CVE) and Power (RMS)
% CVE = 0.523 reflects Gaussian noise
% CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
% CVE > 0.523 reflects phasic activity (e.g. avalanches)
%
% Author: nicolas.gravel@gmail.com
%
% Reference:
% Hidalgo VM, Diaz J, Mpodozis J, Letelier JC. Envelope Analysis of 
% the Human Alpha Rhythm Reveals EEG Gaussianity. IEEE Trans Biomed Eng. 
% 2023 Apr;70(4):1242-1251. 
% doi: 10.1109/TBME.2022.3213840. Epub 2023 Mar 21. PMID: 36223351.

clear all
close all


conditions = ['A', 'B'];
referencing = 1;
switch referencing
    case 1
        load('test_results.mat')
    case 2
        load('test_results_ref.mat')
end

%
tic
for freq_i = 1:3
    
    
    %% Plot RMS envelope as function of CVE
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2)-400*freq_i+200 800, 400]); % Set plot size
    set(gcf, 'color', 'w'); % Set figure background
    
    for cond_i = 1:2
        subplot(1,2,cond_i);
        
        
        X =  results.CVE{freq_i, cond_i}(:); % Coefficient of variation of the envelope
        Y =  results.POW{freq_i, cond_i}(:); % Envelope amplitude
        
        minP  = min(Y(:));
        maxP  = max(Y(:));
        [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/25:maxP});
        [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
        c.LineWidth = 0.0001;
        %pcolor(bin_centers{1},bin_centers{2}, counts');
        
        line([0.523 0.523],[minP maxP],'Color',[1 0 0])
        
        hold on
        axis square;
        colorbar;
        
        %maxCount = 1000; %400;
        %caxis([0 maxCount]);
        %shading interp;
        title(conditions(cond_i),'FontSize', 14);
        %title('baseline time-course','FontSize', 16);
        xlabel('CVE','FontSize', 14);
        ylabel('RMS envelope','FontSize', 14);
        set(gca, 'FontSize', 14);
        set(gca,'LineWidth',1)
        
        clear X Y
    end
    colormap(flipud(bone))
    set(gcf, 'color', 'w'); % Set figure background
    folder = pwd
    switch referencing
        case 1
            fname = [folder '/rmsenv-cve_noref_freq_'  num2str(freq_i) '.png'];
        case 2
            fname = [folder '/rmsenv-cve_ref_freq_'  num2str(freq_i) '.png'];
    end
    %print(gcf, fname, '-dpng', '-r150', '-painters')
    
    %% Plot phase coherence as function of CVE
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2)-400*freq_i+200 800, 400]); % Set plot size
    set(gcf, 'color', 'w'); % Set figure background
    
    for cond_i = 1:2
        subplot(1,2,cond_i);

        
        X =  results.CVE{freq_i, cond_i}(:); % Coefficient of variation of the envelope
        Y =  results.MTST_ch{freq_i, cond_i}(:); % Envelope amplitude
        
        minP  = min(Y(:));
        maxP  = max(Y(:));
        [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/25:maxP});
        [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
        c.LineWidth = 0.0001;
        %pcolor(bin_centers{1},bin_centers{2}, counts');
        
        line([0.523 0.523],[minP maxP],'Color',[1 0 0])
        
        hold on
        axis square;
        colorbar;
        
        %maxCount = 1000; %400;
        %caxis([0 maxCount]);
        %shading interp;
        title(conditions(cond_i),'FontSize', 14);
        %title('baseline time-course','FontSize', 16);
        xlabel('CVE','FontSize', 14);
        ylabel(['{\phi}' ' coherence (K)'],'FontSize',14);
        set(gca, 'FontSize', 14);
        set(gca,'LineWidth',1)
        
        clear X Y
    end
    colormap(flipud(bone))
    set(gcf, 'color', 'w'); % Set figure background
    folder = pwd
    switch referencing
        case 1
            fname = [folder '/phicoh-cve_noref_freq_'  num2str(freq_i) '.png'];
        case 2
            fname = [folder '/phicoh-cve_ref_freq_'  num2str(freq_i) '.png'];
    end
    %print(gcf, fname, '-dpng', '-r150', '-painters')
    
    
    %% Plot phase coherence as function of RMS envelope
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2)-400*freq_i+200 800, 400]); % Set plot size
    set(gcf, 'color', 'w'); % Set figure background
    
    for cond_i = 1:2
        subplot(1,2,cond_i);
        
        
        X =  results.POW{freq_i, cond_i}(:); % Coefficient of variation of the envelope
        Y =  results.MTST_ch{freq_i, cond_i}(:); % Envelope amplitude
        minX  = min(X(:));
        maxX  = max(X(:));
        minY  = min(Y(:));
        maxY  = max(Y(:));
        [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/25:maxP});
        [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
        c.LineWidth = 0.0001;
        %pcolor(bin_centers{1},bin_centers{2}, counts');
        
        line([minX maxX],[minY maxY],'Color',[1 0 0])
        
        hold on
        axis square;
        colorbar;
        
        %maxCount = 1000; %400;
        %caxis([0 maxCount]);
        %shading interp;
        title(conditions(cond_i),'FontSize', 14);
        %title('baseline time-course','FontSize', 16);
        xlabel('RMS envelope','FontSize', 14);
        ylabel(['{\phi}' ' coherence (K)'],'FontSize',14);
        set(gca, 'FontSize', 14);
        set(gca,'LineWidth',1)
        
        clear X Y
    end
    colormap(flipud(bone))
    set(gcf, 'color', 'w'); % Set figure background
    folder = pwd
    switch referencing
        case 1
            fname = [folder '/phicoh-rmsenv_noref_freq_'  num2str(freq_i) '.png'];
        case 2
            fname = [folder '/phicoh-rmsenv_ref_freq_'  num2str(freq_i) '.png'];
    end
    %print(gcf, fname, '-dpng', '-r150', '-painters')
    
end
toc


