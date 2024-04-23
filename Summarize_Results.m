%% Using the Coefficient of Variation (CVE) and RMS envelope to assess Rhythm Gaussianity
% CVE = 0.523 reflects Gaussian noise
% CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
% CVE > 0.523 reflects phasic activity (e.g. avalanches)
%
% Authors:
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


clear; close; clc; %close all

%load('RESULT.mat')



%% Load results
for condition = 1 : length(list.conditions)
    fprintf('\n Printing %s:\n', list.conditions{condition})
    data = eval(['results.' list.conditions{condition}]);
    for file = 1: length(data)
        CVE{condition,file}  =  data{file}.cve(:); % CVE
        RMS{condition,file}  =  data{file}.pow(:); % RMS envelope
        COH{condition,file} =  data{file}.kurVar_ch(:); 
        KUR{condition,file}  =  data{file}.kur_all(:); 
    end
end

%% Compute differences in RMS and CVE
% RMS
Y_a =  vertcat(RMS{1,:}); Y_b =  vertcat(RMS{2,:});
[H0, pval] = ttest2(Y_a,Y_b,'tail','left','alpha',0.001)  % RMS differences are significative
% CVE
X_a =  vertcat(CVE{1,:}); X_b =  vertcat(CVE{2,:});    
[H0,pval] = ttest2(X_a,X_b,'tail','left','alpha',0.001)  % CVE  differences are not significative

%% Compute shifts in RMS and CVE for succesive temporal windows
for i = 1:length(Y_a)-1; dY_a(i) = Y_a(i) - Y_a(i+1); end
for i = 1:length(Y_b)-1; dY_b(i) = Y_b(i) -Y_b(i+1); end
mean(dY_a)
mean(dY_b)
for i = 1:length(X_a)-1; dX_a(i) = X_a(i) - X_a(i+1); end
for i = 1:length(X_b)-1; dX_b(i) = X_b(i) -X_b(i+1); end
mean(dX_a)
mean(dX_b)

%% Compute sequences of RMS shift sign: decreases or increases in power
excursion_length = 3 % sequence length as number of succesive windows, here 9 seconds
%% Go ahead...
% Conditions A
Ex = dY_a < 0;
max_excursionLength = max(diff(find([1,diff(Ex),1]))); % empirical 
runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
%Excursions = Ex(starts(find(runs == excursion_length)));
Excursions = Ex(starts(runs > excursion_length));
decreases_A   = length(nonzeros(Excursions));
increases_A   = length(nonzeros(~Excursions));
% Percentage of events across all windows with length >= excursion_length
power_dyn_A(1,:) = (100/length(Excursions))*decreases_A; % = 44.9541% of alpha power decreasing events longer than excursion_length
power_dyn_A(2,:) = (100/length(Excursions))*increases_A;  % = 55.0459% of alpha power increasing events longer than excursion_length
%(100/length(Excursions))*decreases_A  + (100/length(Excursions))*increases_A % Sanity check 
power_dyn_A

% Conditions B
Ex = dY_b < 0;
max_excursionLength = max(diff(find([1,diff(Ex),1]))); % empirical 
runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
%Excursions = Ex(starts(find(runs == excursion_length)));
Excursions = Ex(starts(runs > excursion_length));
decreases_B   = length(nonzeros(Excursions));
increases_B   = length(nonzeros(~Excursions));
% Percentage of events across all windows with length >= excursion_length
power_dyn_B(1,:) = (100/length(Excursions))*decreases_B; % = 40.5172% of alpha power decreasing events longer than excursion_length
power_dyn_B(2,:) = (100/length(Excursions))*increases_B;  % = 59.4829% of alpha power increasing events longer than excursion_length
%(100/length(Excursions))*decreases_B  + (100/length(Excursions))*increases_B % Sanity check 
power_dyn_B

%% Compute sequences of CVE sign shift: departures from rhythm Gaussianity
%% Go ahead...
% Conditions A
Ex = dX_a > 0;
max_excursionLength = max(diff(find([1,diff(Ex),1])));
runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
%Excursions = Ex(starts(find(runs == excursion_length)));
Excursions = Ex(starts(runs > excursion_length));
towards_rhythm_A  = length(nonzeros(Excursions));
towards_pulses_A   = length(nonzeros(~Excursions));
% Percentage of events across all windows with length > excursion_length
cve_dyn_A(1,:) =  (100/length(Excursions))*towards_rhythm_A; % = 52.1739% of "towards rhythm events" longer than excursion_length
cve_dyn_A(2,:) =  (100/length(Excursions))*towards_pulses_A; % = 47.8261% of "towards pulses" events longer than excursion_length
%(100/length(Excursions))*towards_rhythm_A  + (100/length(Excursions))*towards_pulses_A % Sanity check 
cve_dyn_A

% Conditions B
Ex = dX_b > 0;
max_excursionLength = max(diff(find([1,diff(Ex),1])));
runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
%Excursions = Ex(starts(find(runs == excursion_length)));
Excursions = Ex(starts(runs > excursion_length));
towards_rhythm_B  = length(nonzeros(Excursions));
towards_pulses_B   = length(nonzeros(~Excursions));
% Percentage of events across all windows with length > excursion_length
cve_dyn_B(1,:) =  (100/length(Excursions))*towards_rhythm_B; % = 48.3871% of "towards rhythm events" longer than excursion_length
cve_dyn_B(2,:) =  (100/length(Excursions))*towards_pulses_B; % = 51.6129% of "towards pulses" events longer than excursion_length
%(100/length(Excursions))*towards_rhythm_B  + (100/length(Excursions))*towards_pulses_B % Sanity check 
cve_dyn_B

%% Occupancy: windows below or above 5% and 95% percentiles in CVE for each conditions

Gaussian_CVE = sqrt((4-pi)/pi) % Gaussian CVE
cut_left   = prctile([dX_a,dX_b],5); cut_right = prctile([dX_a,dX_b],95); % percentile cut-offs
thr_left   = Gaussian_CVE + cut_left; thr_right = Gaussian_CVE + cut_right; % thresholds
rhythmicity_a = X_a(X_a <= thr_left); rhythmicity_b = X_b(X_b <= thr_left); % rhythmic windows
phasic_a = X_a(X_a >= thr_right); phasic_b = X_b(X_b >= thr_right); % phasic windows
rhythmic_windows(1,:)= length(rhythmicity_a); rhythmic_windows(2,:)= length(rhythmicity_b);
rhythmic_occupancy(1,:) = (100/length(X_a))*rhythmic_windows(1,:);
rhythmic_occupancy(2,:) = (100/length(X_b))*rhythmic_windows(2,:);
rhythmic_occupancy % percentage of rhytmic windows for each condition
phasic_windows(1,:) = length(phasic_a); phasic_windows(2,:) = length(phasic_b);
phasic_occupancy(1,:) = (100/length(X_a))*phasic_windows(1,:);
phasic_occupancy(2,:) = (100/length(X_b))*phasic_windows(2,:);
phasic_occupancy % percentage of phasic windows for each condition

% Gaussian Noise
Gaussian_noise = wgn(100000,2,0);
window             = 1000; overlap = 500;  dts     = 100;  step     = window - overlap;
win_i = 1;
swin = (win_i*step/dts+1):(win_i*step/dts+window/dts+1)
dim                   = size(Gaussian_noise); tpoints = dim(1);
W                      = round((tpoints-window/dts-dts)/(window-overlap)*dts);
W
f_sampling        = 250; f_low  = 8; f_high = 13; % Hz
f_nrm_low  = f_low /(f_sampling/2); f_nrm_high  = f_high /(f_sampling/2);
% Determine filter coefficients:
[z,p,k] = butter(4,[f_nrm_low f_nrm_high],'bandpass');
% Convert to zero-pole-gain filter parameter (recommended)
sos        = zp2sos(z,p,k);
sig_flt    = sosfilt(sos,Gaussian_noise); size(sig_flt); % apply filter
hbert     = hilbert(sig_flt); size(hbert);
envelope        = abs(real(hbert))';
RMS_envelope = zeros(dim(2),W);
CVE                 = zeros(dim(2),W);
for win_i = 0:W-1
    for ch_i = 1:dim(2)        
        X = envelope(ch_i,(win_i*step/dts+1):(win_i*step/dts+window/dts+1));
        RMS_envelope(ch_i, win_i+1) = rms(X);
        CVE(ch_i, win_i+1) = std(X)/mean(X);
    end
end

%% CVE of Gaussian noise
mean(CVE(:))
Gaussian_CVE = sqrt((4-pi)/pi) % Gaussian CVE
thr_left   = prctile(CVE(:),5)
thr_right = prctile(CVE(:),95) 
%thr_left   = Gaussian_CVE-prctile(CVE(:),1)
%thr_right = Gaussian_CVE+prctile(CVE(:),1)
rhythmicity_a = X_a(X_a <= thr_left); rhythmicity_b = X_b(X_b <= thr_left); % rhythmic windows
phasic_a = X_a(X_a >= thr_right); phasic_b = X_b(X_b >= thr_right); % phasic windows
rhythmic_windows(1,:)= length(rhythmicity_a); rhythmic_windows(2,:)= length(rhythmicity_b);
rhythmic_occupancy(1,:) = (100/length(X_a))*rhythmic_windows(1,:);
rhythmic_occupancy(2,:) = (100/length(X_b))*rhythmic_windows(2,:);
rhythmic_occupancy % percentage of rhytmic windows for each condition
phasic_windows(1,:) = length(phasic_a); phasic_windows(2,:) = length(phasic_b);
phasic_occupancy(1,:) = (100/length(X_a))*phasic_windows(1,:);
phasic_occupancy(2,:) = (100/length(X_b))*phasic_windows(2,:);
phasic_occupancy % percentage of phasic windows for each condition




%% Metastability
kuramoto = KUR(1,:);
metastability(1,:) = std(reshape(cell2mat(vertcat(kuramoto{:})).',1,[]));
kuramoto = KUR(2,:);
metastability(2,:) = std(reshape(cell2mat(vertcat(kuramoto{:})).',1,[]));
metastability


%% CVE histograms
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-600 pos(2) 600, 250]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  Y_a; %results{1}.cve(:);
Xb =  Y_b; %results{2}.cve(:);
C=min([Xa;Xb]):(max([Xa;Xb])/100):max([Xa;Xb]);
[Ns]=histc(Xa,C);
relativefreqNs = Ns ./ sum(Ns);
b(1) = bar(C,relativefreqNs,'histc');
hold on
[Ns]=histc(Xb,C);
relativefreqNs = Ns ./ sum(Ns);
b(2) = bar(C,relativefreqNs,'histc');
xlabel('RMS','FontSize', 14);
ylabel('Relative frequency','FontSize', 14);
%xlim([0.5 1]);
b = findobj(gca,'Type','patch');
set(b(1),'FaceColor', 'k','EdgeColor', 'k','facealpha',0.5,'edgealpha',0);
set(b(2),'FaceColor', 'r','EdgeColor', 'r','facealpha',0.5,'edgealpha',0);
set(gca, 'FontSize', 14);
set(gca,'LineWidth',1)
folder = pwd;
fname = [folder '/hist_rms_alpha_all.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')

figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 600, 250]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  X_a; %results{1}.cve(:);
Xb =  X_b; %results{2}.cve(:);
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
xlim([0 1]);
b = findobj(gca,'Type','patch');
set(b(1),'FaceColor', 'k','EdgeColor', 'k','facealpha',0.5,'edgealpha',0);
set(b(2),'FaceColor', 'r','EdgeColor', 'r','facealpha',0.5,'edgealpha',0);
set(gca, 'FontSize', 14);
set(gca,'LineWidth',1)
folder = pwd;
fname = [folder '/hist_cve_alpha_all.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')

% CVE shifts
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)+600 pos(2) 600, 250]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  dX_a';
Xb =  dX_b';
C=min([Xa;Xb]):(max([Xa;Xb])/100):max([Xa;Xb]);
[Ns]=histc(Xa,C);
relativefreqNs = Ns ./ sum(Ns);
b(1) = bar(C,relativefreqNs,'histc');
hold on
[Ns]=histc(Xb,C);
relativefreqNs = Ns ./ sum(Ns);
b(2) = bar(C,relativefreqNs,'histc');
xlabel('CVE shifts','FontSize', 14);
ylabel('Relative frequency','FontSize', 14);
%xlim([0.5 1]);
b = findobj(gca,'Type','patch');
set(b(1),'FaceColor', 'k','EdgeColor', 'k','facealpha',0.5,'edgealpha',0);
set(b(2),'FaceColor', 'r','EdgeColor', 'r','facealpha',0.5,'edgealpha',0);
set(gca, 'FontSize', 14);
set(gca,'LineWidth',1)
folder = pwd;
fname = [folder '/hist_cve_shifts_all.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')

% Kuramoto variance histograms
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-400 pos(2)-100 600, 250]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  vertcat(COH{1,:}); % Kuramoto
Xb =  vertcat(COH{2,:}); % Kuramoto
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
fname = [folder '/kurvar_alpha_all.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')


%% Plot RMS envelope as function of CVE
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-400 pos(2)-100 800, 400]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
% The set both axis (conditions) with the same range
X_1 =  vertcat(CVE{1,:}); % CVE
Y_1 =  vertcat(RMS{1,:}); % RMS envelope
X_2 =  vertcat(CVE{2,:}); % CVE
Y_2 =  vertcat(RMS{2,:});  % RMS envelope
minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
binfactor = 25;
maxCount = 250;
for cond_i = 1:2
    subplot(1,2,cond_i);
    X =  vertcat(CVE{cond_i,:}); % CVE
    Y =  vertcat(RMS{cond_i,:}); % RMS envelope
    minP  = minY; %min(Y(:));
    maxP  = maxY; %max(Y(:));
    bin = length(Y)/binfactor;
    [counts, bin_centers] = hist3([X,Y],'Ctrs',{0:0.05:1 0:maxP/binfactor:maxP});
    [~,c] = contourf(bin_centers{1},bin_centers{2}, counts');
    c.LineWidth = 0.0001;
    %pcolor(bin_centers{1},bin_centers{2}, counts');
    line([0.523 0.523],[minY maxY],'Color',[1 0 0])
    hold on
    xlim([0 1]);
    %ylim([minY maxY]);
    ylim([minY 2]);
    axis square;
    colorbar;
    caxis([0 maxCount]);
    %title([list.conditions(cond_i) ' ( #' num2str(length(results{cond_i}.time)) ')'],'FontSize', 14);
    xlabel('CVE','FontSize', 14);
    ylabel('RMS envelope','FontSize', 14);
    set(gca,'FontSize', 14);
    set(gca,'LineWidth',1)
    clear X Y
end
colormap(flipud(bone))
folder = pwd;
fname = [folder '/rmsenv-cve_alpha_all.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')


%% Plot phase coherence as function of RMS envelope
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)+400 pos(2)+200 800, 400]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
X_1 =  vertcat(COH{1,:}); % Phase coherence
Y_1 =  vertcat(RMS{1,:}); % RMS envelope
X_2 =  vertcat(COH{2,:}); % Phase coherence
Y_2 =  vertcat(RMS{2,:});  % RMS envelope
minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
binfactor = 25;
for cond_i = 1:2
    subplot(1,2,cond_i);
    %X =  results{cond_i}.kurVar_ch(:);  % Phase coherence
    %Y =  results{cond_i}.pow(:);        % RMS envelope
    X =  vertcat(COH{cond_i,:}); % Phase coherence
    Y =  vertcat(RMS{cond_i,:});  % RMS envelope
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
    %title([conditions(cond_i) ' ( #' num2str(length(results{cond_i}.time)) ')'],'FontSize', 14);
    ylabel('RMS envelope','FontSize', 14);
    xlabel(['std({\Delta}{\phi} coherence (R) )'],'FontSize', 14);
    set(gca, 'FontSize', 14);
    set(gca,'LineWidth',1)
    xlim([0 maxX]);
    ylim([0 maxY]);
    clear X Y
end
colormap(flipud(bone))
folder = pwd
fname = [folder '/kurvar-rms_alpha_all.png'];
print(gcf, fname, '-dpng', '-r150', '-painters')


 