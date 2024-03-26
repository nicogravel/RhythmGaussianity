%% Compute Coefficient of Variation (CVE) and Power (RMS)
% CVE = 0.523 reflects Gaussian noise
% CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
% CVE > 0.523 reflects phasic activity (e.g. avalanches)
% (Koiti Motokawa, 1943)

%% Band-pass filter
freqs_low = [2, 10, 30]; %2 % Hz
freqs_high = [8, 20, 60]; %12; % Hz   (Hidalgo et al. use 45 Hz)
conditions = ['A', 'B'];

referencing = 1
switch referencing
    case 1
        data      = data_unref
    case 2
        data      = data_ref
end

for freq_i = 1:3
    
    
    %% Power Dynamics
    f_low  = freqs_low(freq_i);
    f_high = freqs_high(freq_i);
    
    f_sampling = 250; % 1kHz
    f_nrm_low   = f_low /(f_sampling/2);
    f_nrm_high  = f_high /(f_sampling/2);
    % determine filter coefficients:
    [z,p,k] = butter(4,[f_nrm_low f_nrm_high],'bandpass');
    % convert to zero-pole-gain filter parameter (recommended)
    sos = zp2sos(z,p,k);
    % determine sliding window
    window  = 6000;
    overlap = 2000;
    dts     = 500;
    step    = window - overlap;
    maxP = 1 %max(Y(:));
    maxCount = 400;
    trajectory = 0
    
    % plot Power as function of CVE
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2)-400*freq_i+400 800, 400]); % Set plot size
    set(gcf, 'color', 'w'); % Set figure background
    
    for cond_i = 1:2
        subplot(1,2,cond_i);
        
        
        % Compute CVE
        %figure,
        sig = zscore(data{cond_i}(:)); % load and zscore data
        sig_flt = sosfilt(sos,sig); % apply filter
        tpoints = size(sig_flt,1);
        % determine sliding window
        %window  = 4000;
        %overlap = 2000;
        %dts     = 20;
        step    = window - overlap;
        W  = round((tpoints-window/dts-dts)/(window-overlap)*dts);
        W;
        
        
        % apply Hilbert transform
        h = hilbert(sig_flt);
        envelope = abs(real(h)); % envelope
        size(envelope)
        
        % phase differences
        %sig_ch = zscore(data{cond_i}');
        % apply filter
        %sig_flt = sosfilt(sos,sig_ch);
        %h = hilbert(sig_flt);
        phase    = angle(h); % angle
        size(phase)
        
        for win_i = 0:W-1
            
            % envelope
            wenv = envelope((win_i*step/dts+1):(win_i*step/dts+window/dts+1));
            sd = std(wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10))));
            me = mean(wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10))));
            cve(:,win_i+1) = sd/me;
            pow(:,win_i+1) = rms(wenv);
            
            
        end
        
        results.CVE{cond_i} = cve;
        results.POW{cond_i} = pow;
        
        clear cve pow envelope wenv
        % Histogram of power
        %figure,hist(cve(:),100)
        %figure,hist(pow(:),100)
        
        X =  results.CVE{cond_i}(:); % Coefficient of variation of the envelope
        Y =  results.POW{cond_i}(:); % Envelope amplitude
        %maxP = 1;
        maxP  = max(Y(:));
        [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/25:maxP});
        [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
        c.LineWidth = 0.0001;
        %pcolor(bin_centers{1},bin_centers{2}, counts');
        
        line([0.523 0.523],[0 maxP],'Color',[1 0 0])
        
        hold on
        
        
        % For plotting the numbers and eventually derive trajectory
        if trajectory == 1
            for win_i = 0:W-1
                x = results.CVE{cond_i}(:,win_i+1);
                y =  results.POW{cond_i}(:,win_i+1);
                text(x,y,num2str(win_i),'FontSize',14);
                hold on
            end
        end
        
        
        axis square;
        colorbar;
        
        maxCount = 1000; %400;
        %caxis([0 maxCount]);
        %shading interp;
        title(conditions(cond_i),'FontSize', 14);
        %title('baseline time-course','FontSize', 16);
        xlabel('CVE','FontSize', 14);
        ylabel('RMS envelope','FontSize', 14);
        set(gca, 'FontSize', 14);
        set(gca,'LineWidth',1)
        
    end
    colormap(flipud(bone))
    set(gcf, 'color', 'w'); % Set figure background
    folder = pwd
    switch referencing
        case 1
            fname = [folder '/rmsenvDyn_noref_freq_'  num2str(freq_i) '.png'];
        case 2
            fname = [folder '/rmsenvDyn_Laplace_freq_'  num2str(freq_i) '.png'];
    end
    print(gcf, fname, '-dpng', '-r150', '-painters')
    
    
    
    % plot standard deviation of phase as a function of CVE
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)+800 pos(2)-400*freq_i+400 800, 400]); % Set plot size
    set(gcf, 'color', 'w'); % Set figure background
    
    for cond_i = 1:2
        subplot(1,2,cond_i);
        
        
        % Compute CVE
        %figure,
        sig = zscore(data{cond_i}(:)); % load and zscore data
        sig_flt = sosfilt(sos,sig); % apply filter
        tpoints = size(sig_flt,1);
        % determine sliding window
        %window  = 4000;
        %overlap = 2000;
        %dts     = 20;
        step    = window - overlap;
        W  = round((tpoints-window/dts-dts)/(window-overlap)*dts);
        W;
        
        
        % apply Hilbert transform
        h = hilbert(sig_flt);
        envelope = abs(real(h)); % envelope
        size(envelope)
        
        % phase differences
        %sig_ch = zscore(data{cond_i}');
        % apply filter
        %sig_flt = sosfilt(sos,sig_ch);
        %h = hilbert(sig_flt);
        phase    = angle(h); % angle
        size(phase)
        
        
        for win_i = 0:W-1
            
            % envelope
            wenv = envelope((win_i*step/dts+1):(win_i*step/dts+window/dts+1));
            sd = std(wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10))));
            me = mean(wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10))));
            cve(:,win_i+1) = sd/me;
            pow(:,win_i+1) = rms(wenv);
            
            % phase
            wphase = phase((win_i*step/dts+1):(win_i*step/dts+window/dts+1));
            phi    = wphase(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10)));
            dPhi = bsxfun(@minus,phi,phi');
            D = wrapToPi(dPhi(find(tril(ones(size(dPhi)),-1))));
            kurDev(:,win_i+1) = std(abs(real(exp(1i*D))));
            
            %clear linkOrder meta dPhi D
            %for t = 1:size(phi,1) % loop over time
            %    dPhi = bsxfun(@minus,phi(t,:),phi(t,:)');
            %    D = wrapToPi(dPhi(find(tril(ones(size(deltaPh)),-1))));
            %    %abs(real(sum(exp(1i*phi(t,:))))/size(dPhi,1));
            %    linkOrder(t) = abs(sum(exp(1i*phi(t,:))))/size(dPhi,1);
            %end
            %meta(:,win_i+1) = std(linkOrder);
            
        end
        
        results.CVE{cond_i} = cve;
        results.MET{cond_i} = kurDev;
        
        clear cve pow envelope wenv kurDev
        % Histogram of power
        %figure,hist(cve(:),100)
        %figure,hist(pow(:),100)
        
        X =  results.CVE{cond_i}(:); % Coefficient of variation of the envelope
        Y =  results.MET{cond_i}(:); % Envelope amplitude
        minP  = min(Y(:));
        maxP  = max(Y(:));
        [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 minP:maxP/25:maxP});
        [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
        c.LineWidth = 0.0001;
        %pcolor(bin_centers{1},bin_centers{2}, counts');
        
        line([0.523 0.523],[minP maxP],'Color',[1 0 0])
        
        
        axis square;
        colorbar;
        
        maxCount = 1000; %400;
        %caxis([0 maxCount]);
        %shading interp;
        title(conditions(cond_i),'FontSize', 14);
        %title('baseline time-course','FontSize', 16);
        xlabel('CVE','FontSize', 14);
        ylabel(['{\phi}' ' coherence (K)'],'FontSize',14);
        set(gca, 'FontSize', 14);
        set(gca,'LineWidth',1)
        
    end
    colormap(flipud(bone))
    set(gcf, 'color', 'w'); % Set figure background
    folder = pwd
    switch referencing
        case 1
            fname = [folder '/phaseDyn_noref_freq_'  num2str(freq_i) '.png'];
        case 2
            fname = [folder '/phaseDyn_Laplace_freq_'  num2str(freq_i) '.png'];
    end
    print(gcf, fname, '-dpng', '-r150', '-painters')
    
    
    
    % plot standard deviation of phase as a function of power
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)+800*2 pos(2)-400*freq_i+400 800, 400]); % Set plot size
    set(gcf, 'color', 'w'); % Set figure background
    
    for cond_i = 1:2
        subplot(1,2,cond_i);
        
        
        % Compute CVE
        %figure,
        sig = zscore(data{cond_i}(:)); % load and zscore data
        sig_flt = sosfilt(sos,sig); % apply filter
        tpoints = size(sig_flt,1);
        % determine sliding window
        %window  = 4000;
        %overlap = 2000;
        %dts     = 20;
        step    = window - overlap;
        W  = round((tpoints-window/dts-dts)/(window-overlap)*dts);
        W;
        
        
        % apply Hilbert transform
        h = hilbert(sig_flt);
        envelope = abs(real(h)); % envelope
        size(envelope)
        
        % phase differences
        %sig_ch = zscore(data{cond_i}');
        % apply filter
        %sig_flt = sosfilt(sos,sig_ch);
        %h = hilbert(sig_flt);
        phase    = angle(h); % angle
        size(phase)
        
        for win_i = 0:W-1
            
            % envelope
            wenv = envelope((win_i*step/dts+1):(win_i*step/dts+window/dts+1));
            sd = std(wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10))));
            me = mean(wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10))));
            cve(:,win_i+1) = sd/me;
            pow(:,win_i+1) = rms(wenv);
            
            % phase
            wphase = phase((win_i*step/dts+1):(win_i*step/dts+window/dts+1));
            phi    = wphase(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10)));
            dPhi = bsxfun(@minus,phi,phi');
            D = wrapToPi(dPhi(find(tril(ones(size(dPhi)),-1))));
            kurDev(:,win_i+1) = std(abs(real(exp(1i*D))));
            
            %clear linkOrder meta dPhi D
            %for t = 1:size(phi,1) % loop over time
            %    dPhi = bsxfun(@minus,phi(t,:),phi(t,:)');
            %    D = wrapToPi(dPhi(find(tril(ones(size(deltaPh)),-1))));
            %    %abs(real(sum(exp(1i*phi(t,:))))/size(dPhi,1));
            %    linkOrder(t) = abs(sum(exp(1i*phi(t,:))))/size(dPhi,1);
            %end
            %meta(:,win_i+1) = std(linkOrder);
            
        end
        
        results.CVE{cond_i} = cve;
        results.MET{cond_i} = kurDev;
        
        clear cve pow envelope wenv kurDev
        % Histogram of power
        %figure,hist(cve(:),100)
        %figure,hist(pow(:),100)
        
        X =  results.POW{cond_i}(:); % Coefficient of variation of the envelope
        Y =  results.MET{cond_i}(:); % Envelope amplitude
        minX  = min(X(:));
        maxX  = max(X(:));
        minY  = min(Y(:));
        maxY  = max(Y(:));
        [counts, bin_centers] = hist3([X,Y],'Ctrs',{minX:maxX/25:maxX minP:maxP/25:maxP});
        [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
        c.LineWidth = 0.0001;
        %pcolor(bin_centers{1},bin_centers{2}, counts');
        
        line([minX maxX],[minY maxY],'Color',[1 0 0])
        line([minX maxX],[minY maxY],'Color',[1 0 0])
        
        
        axis square;
        colorbar;
        
        maxCount = 1000; %400;
        %caxis([0 maxCount]);
        %shading interp;
        title(conditions(cond_i),'FontSize', 14);
        %title('baseline time-course','FontSize', 16);
        xlabel('RMS envelope','FontSize', 14);
        ylabel(['{\phi}' ' coherence (K)'],'FontSize',14);
        set(gca, 'FontSize', 14);
        set(gca,'LineWidth',1)
        
    end
    colormap(flipud(bone))
    set(gcf, 'color', 'w'); % Set figure background
    folder = pwd
    switch referencing
        case 1
            fname = [folder '/ph-rmsnenvDyn_noref_freq_'  num2str(freq_i) '.png'];
        case 2
            fname = [folder '/ph-rmsenvDyn_Laplace_freq_'  num2str(freq_i) '.png'];
    end
    print(gcf, fname, '-dpng', '-r150', '-painters')
    
    
    
    
end


