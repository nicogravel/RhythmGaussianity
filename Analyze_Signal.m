%% Using the Coefficient of Variation (CVE) and RMS envelope to assess Rhythm Gaussianity
% CVE = 0.523 reflects Gaussian noise
% CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
% CVE > 0.523 reflects phasic activity (e.g. avalanches)
%
% Author:
% Nicolas Gravel <nicolas.gravel@gmail.com>
%
% References:
%
% Hidalgo VM, Diaz J, Mpodozis J, Letelier JC. Envelope Analysis of
% the Human Alpha Rhythm Reveals EEG Gaussianity. IEEE Trans Biomed Eng.
% 2023 Apr;70(4):1242-1251.
% doi: 10.1109/TBME.2022.3213840. Epub 2023 Mar 21. PMID: 36223351.
%
% Hidalgo VM, Letelier JC, Diaz J. The amplitude modulation pattern of 
% Gaussian noise is a fingerprint of Gaussianity.  	
% https://doi.org/10.48550/arXiv.2203.16253


clear; clc; %close all

addpath /Users/Nicolas/Downloads/fieldtrip-20231025
addpath /Users/Nicolas/Documents/MATLAB/CFmrVista/vistasoft/cc-pRF/MISC/CircStats
ft_defaults()


%% Analyse data
% These are the two datasets from CZECH institute
files = {'sampleAt1_5hA.edf', 'sampleAt1_5hB.edf'}
conditions        = ['A', 'B'];
fband_of_interest = [8 13];
phasediff_histos  = 0;
win_length        = 24;  % in seconds
overlap           = 0.9; % fraction

for cond_i = 1:2
    
    % Load data
    cfg              = [];
    cfg.dataset      = files{cond_i};
    cfg.channel      = 'EEG'; % this is the default
    data_eeg         = ft_preprocessing(cfg);
    data_eeg.fsample = 1000;
    
    % Re-reference
    cfg            = [];
    cfg.reref      = 'yes';
    cfg.refmethod  = 'avg';
    cfg.refchannel = 'all';
    data_ref       = ft_preprocessing(cfg, data_eeg)
    clear data_eeg
    
    % Shift and scale channels
    cfg            = [];
    cfg.demean     = 'no';
    %cfg.method     = 'perchannel' 
    cfg.method     = 'acrosschannel'
    cfg.channel    = [72, 173, 114, 119, 168]
    %cfg.channel    = [114, 119, 168]
    %cfg.channel    = [36, 72, 114, 224, 173, 168]
    %cfg.channel    = [244, 36, 72, 114, 234, 224, 173, 168]
    [data_norm]    = ft_channelnormalise(cfg, data_ref)
    clear data_ref
    
    % Filter
    cfg            = [];
    cfg.bpfilter   = 'yes'
    cfg.bpfreq     = fband_of_interest;
    cfg.bpfiltord  = 4;
    cfg.bpfilttype = 'but';
    data_filt      = ft_preprocessing(cfg, data_norm)
    clear data_norm
    
    % Compute Hilbert analytical signal
    cfg            = [];
    cfg.hilbert    = 'complex';
    complex_data   = ft_preprocessing(cfg, data_filt)
    clear data_filt
    
    % Sliding window
    cfg           = [];
    cfg.length    = win_length;   % specifies, in seconds, the length of the required snippets
    cfg.overlap   = overlap; % fraction of overlap between snippets (0 = no overlap)
    [data_win]    = ft_redefinetrial(cfg, complex_data)
    clear complex_data
    
    win_N = length(data_win.trial);
    ch_N  = length(data_win.label);
    kurVar_ch = zeros(ch_N,win_N);
    kurVar_all = zeros(ch_N,win_N);
    cve        = zeros(ch_N,win_N);
    pow        = zeros(ch_N,win_N);
    
    for win_i = 1 : win_N
        
        Y        = data_win.trial{win_i};
        envelope = real(abs(Y));
        phase    = wrapToPi(angle(Y));
        
        % phase synchrony from one-to-all channels
        for t = 1:size(phase,2) % loop over time
            ph_diff{t} = bsxfun(@minus,phase(:,t),phase(:,t)');
        end
        phDiff             = wrapToPi(cat(3,ph_diff{:}));
        kuramoto_ch        = squeeze(abs(real(sum(exp(1i*phDiff))))/size(phDiff ,1));
        kurVar_ch(:,win_i) = std(kuramoto_ch,0,2);
        
        % phase synchrony across all channels (excluding self and mirror)
        if phasediff_histos == 1
            binEdges = -pi:pi/25:pi; % phase bin edges
            P        = zeros(size(binEdges,2)-1,size(phase,2));
        end
        
        for t = 1:size(phase,2) % loop over time
            dPh  = phDiff(:,:,t);
            D{t} = dPh(find(tril(ones(size(dPh)),-1)));
            if phasediff_histos == 1
                [N,binEdges] = histcounts(D{t},binEdges);
                P(:,t) = N; % keep distribution (mostly for illustration/scrutiny)
            end
        end
        clear dPh
        if phasediff_histos == 1
            phistog{win_i}    = P; % Saving the spatiotemporal distribution is heavy!
        end
        dPh               = wrapToPi(cell2mat(D));
        kuramoto          = abs(real(sum(exp(1i*dPh))))/size(dPh,1);     % Def 2
        kurVar_all(win_i) = std(kuramoto);
        kur_all{win_i}    = kuramoto;
        
        % CVE and RMS envelope
        %for ch_i = 1 : ch_N
        %    out.cve(ch_i,win_i)  = std(envelope(ch_i,:))/mean(envelope(ch_i,:));
        %    out.pow(ch_i,win_i)  = rms(envelope(ch_i,:));
        %end
        
        % CVE and RMS envelope
        cve(:,win_i)  = std(envelope,0,2)./mean(envelope,2);
        pow(:,win_i)  = rms(envelope,2);
        
    end
    
    %% Keep results
    out.time        = data_win.time;
    out.channels    = data_win.label;
    out.cve         = cve;
    out.pow         = pow;
    out.kur_all     = kur_all;
    out.kurVar_all  = kurVar_all;
    out.kurVar_ch   = kurVar_ch;
    if phasediff_histos == 1
        out.phistog     = phistog;
    end
    results{cond_i} = out;
    clearvars -except results files conditions fband_of_interest phasediff_histos win_length overlap out binEdges
end



%% Plot RMS envelope as function of CVE
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-400 pos(2)-100 800, 400]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
% The set both axis (conditions) with the same range
X_1 =  results{1}.cve(:); % CVE
Y_1 =  results{1}.pow(:); % RMS envelope
X_2 =  results{2}.cve(:); % CVE
Y_2 =  results{2}.pow(:); % RMS envelope
minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
binfactor = 25;
for cond_i = 1:2
    subplot(1,2,cond_i);
    X =  results{cond_i}.cve(:); % CVE
    Y =  results{cond_i}.pow(:); % RMS envelope
    minP  = min(Y(:));
    maxP  = max(Y(:));
    bin = length(Y)/binfactor;
    [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/binfactor:maxP});
    [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
    c.LineWidth = 0.0001;
    %pcolor(bin_centers{1},bin_centers{2}, counts');
    line([0.523 0.523],[minY maxY],'Color',[1 0 0])
    hold on
    xlim([0 1]);
    ylim([minY maxY]);
    axis square;
    colorbar;
    title([conditions(cond_i) ' ( #' num2str(length(results{cond_i}.time)) ')'],'FontSize', 14);
    xlabel('CVE','FontSize', 14);
    ylabel('RMS envelope','FontSize', 14);
    set(gca,'FontSize', 14);
    set(gca,'LineWidth',1)
    clear X Y
end
colormap(flipud(bone))
folder = pwd
fname = [folder '/rmsenv-cve_alpha_v1.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')

% % Plot RMS envelope as function of CVE
% figure,
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1)-400 pos(2)-100 800, 400]); % Set plot size
% set(gcf, 'color', 'w'); % Set figure background
% binX  = 50;
% binY  = 50;
% % The set both axis (conditions) with the same range
% X_1 =  results{1}.cve(:); % CVE
% Y_1 =  results{1}.pow(:); % RMS envelope
% X_2 =  results{2}.cve(:); % CVE
% Y_2 =  results{2}.pow(:); % RMS envelope
% minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
% minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
% for cond_i = 1:2
%     subplot(1,2,cond_i);
%     X =  results{cond_i}.cve(:); % CVE
%     Y =  results{cond_i}.pow(:); % RMS envelope
%     cMax = 1;
%     xLabel = 'CVE';
%     yLabel = 'RMS envelope';
%     %# bin centers (integers)
%     xbins = floor(min(X)):1/binX:ceil(max(X));
%     ybins = floor(min(Y)):1/binY:ceil(max(Y));
%     xNumBins = numel(xbins); yNumBins = numel(ybins);
%     Xi = round( interp1(xbins, 1:xNumBins, X, 'linear', 'extrap') );
%     Yi = round( interp1(ybins, 1:yNumBins, Y, 'linear', 'extrap') );
%     Xi = max( min(Xi,xNumBins), 1);
%     Yi = max( min(Yi,yNumBins), 1);
%     H = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
%     imagesc(xbins, ybins, H), axis on
%     hold on
%     line([0.523 0.523],[minY maxY],'Color',[1 0 0])
%     data=cat(1,X,Y);
%     data=sortrows(data,1)';
%     %caxis([0 cMax]);
%     cb=colorbar;
%     xlabel(cb,'#', 'FontSize', 14);
%     set(gca,'YDir','normal')
%     set(gca, 'FontSize', 14);
%     set(gca,'LineWidth',1)
%     set(gcf, 'color', 'w');
%     xlabel(xLabel,'FontSize', 14);
%     xlim([0 1]);
%     ylim([minY maxY]);
%     ylabel(yLabel,'FontSize', 14);
%     axis square
% end
% colormap(flipud(bone))
% folder = pwd
% fname = [folder '/rmsenv-cve_alpha_v2.png'];
% print(gcf, fname, '-dpng', '-r150', '-painters')


%% Plot phase coherence as function of RMS envelope
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)+400 pos(2)+200 800, 400]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
% The set both axis (conditions) with the same range
X_1 =  results{1}.kurVar_ch(:);   % Phase coherence
Y_1 =  results{1}.pow(:);         % RMS envelope
X_2 =  results{2}.kurVar_ch(:);   % Phase coherence
Y_2 =  results{2}.pow(:);         % RMS envelope
minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
binfactor = 25;
for cond_i = 1:2
    subplot(1,2,cond_i);
    X =  results{cond_i}.kurVar_ch(:);  % Phase coherence
    Y =  results{cond_i}.pow(:);        % RMS envelope
    minP  = min(Y(:));
    maxP  = max(Y(:));
    bin = length(Y)/binfactor;
    [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/binfactor:maxP});
    [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
    c.LineWidth = 0.0001;
    %pcolor(bin_centers{1},bin_centers{2}, counts');
    hold on
    line([0 maxX],[0 maxY],'Color',[1 0 0])
    axis square;
    colorbar;
    title([conditions(cond_i) ' ( #' num2str(length(results{cond_i}.time)) ')'],'FontSize', 14);
    ylabel('RMS envelope','FontSize', 14);
    xlabel(['{\Delta}{\phi}' ' coherence (R)'],'FontSize',14);
    set(gca, 'FontSize', 14);
    set(gca,'LineWidth',1)
    xlim([0 maxX]);
    ylim([0 maxY]);
    clear X Y
end
colormap(flipud(bone))
folder = pwd
fname = [folder '/kurvar-rms_alpha_v1.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')

% % Plot phase coherence as function of RMS envelope
% figure,
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1)-400 pos(2)-100 800, 400]); % Set plot size
% set(gcf, 'color', 'w'); % Set figure background
% binX  = 100;
% binY  = 50;
% X_1 =  results{1}.kurVar_ch(:);   % Phase coherence
% Y_1 =  results{1}.pow(:);         % RMS envelope
% X_2 =  results{2}.kurVar_ch(:);   % Phase coherence
% Y_2 =  results{2}.pow(:);         % RMS envelope
% minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
% minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
% for cond_i = 1:2
%     subplot(1,2,cond_i);
%     X =  results{cond_i}.kurVar_ch(:); % CVE
%     Y =  results{cond_i}.pow(:); % RMS envelope
%     cMax = 1;
%     xLabel = 'CVE';
%     yLabel = ['{\Delta}{\phi}' ' coherence (R)'];
%     %# bin centers (integers)
%     xbins = floor(min(X)):1/binX:ceil(max(X));
%     ybins = floor(min(Y)):1/binY:ceil(max(Y));
%     xNumBins = numel(xbins); yNumBins = numel(ybins);
%     Xi = round( interp1(xbins, 1:xNumBins, X, 'linear', 'extrap') );
%     Yi = round( interp1(ybins, 1:yNumBins, Y, 'linear', 'extrap') );
%     Xi = max( min(Xi,xNumBins), 1);
%     Yi = max( min(Yi,yNumBins), 1);
%     H = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
%     imagesc(xbins, ybins, H), axis on
%     hold on
%     line([0 maxX],[0 maxY],'Color',[1 0 0])
%     data=cat(1,X,Y);
%     data=sortrows(data,1)';
%     %caxis([0 cMax]);
%     cb=colorbar;
%     xlabel(cb,'#', 'FontSize', 14);
%     set(gca,'YDir','normal')
%     set(gca, 'FontSize', 14);
%     set(gca,'LineWidth',1)
%     set(gcf, 'color', 'w');
%     xlabel(xLabel,'FontSize', 14);
%     xlim([0 maxX]);
%     ylim([0 maxY]);
%     clear X Y
%     ylabel(yLabel,'FontSize', 14);
%     axis square
% end
% colormap(flipud(bone))
% folder = pwd
% fname = [folder '/kurvar-rms_alpha_v2.png'];
% print(gcf, fname, '-dpng', '-r150', '-painters')


%% CVE histograms
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-400 pos(2)-100 400, 200]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  results{1}.cve(:);
Xb =  results{2}.cve(:);
C=min([Xa;Xb]):(max([Xa;Xb])/100):max([Xa;Xb]);
[Ns]=histc(Xa,C);
relativefreqNs = Ns ./ sum(Ns);
b(1) = bar(C,relativefreqNs,'histc');
hold on
[Ns]=histc(Xb,C);
relativefreqNs = Ns ./ sum(Ns);
b(2) = bar(C,relativefreqNs,'histc');
xlabel('CVE','FontSize', 14);
ylabel('Relative frequency','FontSize', 14);
%xlim([0.5 1]);
b = findobj(gca,'Type','patch');
set(b(1),'FaceColor', 'k','EdgeColor', 'k','facealpha',0.5,'edgealpha',0);
set(b(2),'FaceColor', 'r','EdgeColor', 'r','facealpha',0.5,'edgealpha',0);
set(gca, 'FontSize', 14);
set(gca,'LineWidth',1)
folder = pwd
fname = [folder '/hist_cve_alpha.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')

%% Kuramoto variance histograms
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-400 pos(2)-100 400, 200]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  results{1}.kurVar_ch(:);
Xb =  results{2}.kurVar_ch(:);
C=0:(max([Xa;Xb])/50):max([Xa;Xb]);
[Ns]=histc(Xa,C);
relativefreqNs = Ns ./ sum(Ns);
b(1) = bar(C,relativefreqNs,'histc');
hold on
[Ns]=histc(Xb,C);
relativefreqNs = Ns ./ sum(Ns);
b(2) = bar(C,relativefreqNs,'histc');
xlabel(['std({\Delta}{\phi} coherence (R) )'],'FontSize', 14);
ylabel('Relative frequency','FontSize', 14);
%xlim([0.5 1]);
b = findobj(gca,'Type','patch');
set(b(1),'FaceColor', 'k','EdgeColor', 'k','facealpha',0.5,'edgealpha',0);
set(b(2),'FaceColor', 'r','EdgeColor', 'r','facealpha',0.5,'edgealpha',0);
set(gca, 'FontSize', 14);
set(gca,'LineWidth',1)
folder = pwd
fname = [folder '/kurvar_alpha.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')


%% Plot spatiotemporal phase differences distribution
binEdges = -pi:pi/25:pi; % phase bin edges
if phasediff_histos == 1
    for cond_i = 1:2
        figure,
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1)-400 pos(2)-600 1200, 1200]); % Set plot size
        set(gcf, 'color', 'w'); % Set figure background
        max_idx = length(results{cond_i}.time)
        win_idx = [30, 50, 60, 80, 90, 100]; % Plot six example windows
        for win_i = 1:length(win_idx)
            subplot(3,2,win_i);
            P = results{cond_i}.phistog{win_i};
            Y = results{cond_i}.kur_all{win_i};
            pcolor(P); shading flat; colormap jet;
            caxis([0 max(max(P))]); % Plot phase diferences time flow histogram
            hold on
            plot(Y*numel(binEdges)+1,'-w','Linewidth',2); % Plot mean phase coherence
            title(['win idx: ' num2str(win_idx(win_i))],'FontSize',8);
            xlabel('Time point (p.d.u.)','FontSize',8,'FontWeight','Bold');
            set(gca, 'FontSize', 8,'FontWeight','Bold', 'LineWidth', 2); % Set plot font size and linewidth
            set(gca,'YTickLabel',[1/5 2/5 3/5 4/5 1]);
            text(size(P,2)+2,0,'-\pi','FontSize',18,'FontWeight','Bold');
            text(size(P,2)+2,25,'0','FontSize',8);
            text(size(P,2)+2,50,'\pi','FontSize',18,'FontWeight','Bold');
            ylabel(['{\Delta}' '{\phi}' ' coherence (R)'],'FontSize',8);
            set(gcf, 'color', 'w');
            set(gca, 'box', 'off');
        end
        folder = pwd
        fname = [folder '/kurHist_alpha_cond' num2str(conditions(cond_i)) '.PNG'];
        print(gcf, fname, '-dpng', '-r300', '-painters')
    end
end

