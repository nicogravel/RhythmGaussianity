%% Itinerancy in the RMS-CVE space 
% To quantify the itinerancy in the RMS-CVE space, the tendency for departures from Gaussianity, 
% we proceeded as follows. First, for each subject and condition, we counted all sequences of contiguous 
% epochs in which the RMS envelope and the CVE were increasing or decreasing. To do so, we defined a
% minimum sequence length (3 consecutive epochs, equivalent to 9 seconds) and identified the start of
% each sequence. Second, we summed the number of contiguous sequences with positive or negative shifts.
% Third, we computed, for each subject and condition, the percentage of the sequences.
%To assess the departures from rhythm Gaussianity, we created a null distribution of CVE values by applying
% the procedure used to process the EEG data to Gaussian noise. This resulted in a distribution of CVE values
% with mean 0.523 (very close to the theoretical value for Gaussian noise: √(4 − π)/π ≈ 0.523). We then used
% the 10% and 90% percentile values of this null distribution as thresholds to define departures from 
% Gaussianity. This allowed us to compute, for each subject and condition, the percentage of epochs in which
% the data departed from Gaussianity.
%
% Author:
% Nicolas Gravel <nicolas.gravel@gmail.com>


clear all; close all; clc;

load('RESULT_withinchan.mat')

excursion_length = 2; % sequence length as number of succesive windows

%threshold = 'Gaussian_Noise'  % just for testing!
threshold = 'Bandpass_Percentiles'

%% Load results
for condition = 1 : length(list.conditions)
    fprintf('\n Printing %s:\n', list.conditions{condition})
    data = eval(['results.' list.conditions{condition}]);
    for file = 1: length(data)
        CVE{condition,file}  =  data{file}.cve(:);  % CVE
        RMS{condition,file}  =  data{file}.pow(:); % RMS envelope
    end
end

%% Compute differences in RMS and CVE
% RMS
Y_a =  vertcat(RMS{1,:}); Y_b =  vertcat(RMS{2,:});
[H0, pval] = ttest2(Y_a,Y_b,'tail','left','alpha',0.001)  % RMS differences are significative
% CVE
X_a =  vertcat(CVE{1,:}); X_b =  vertcat(CVE{2,:});
[H0,pval] = ttest2(X_a,X_b,'tail','left','alpha',0.001)  % CVE  differences are not significative

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


%% Plot bivariate histograms for RMS-CVE
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



%% Compute shifts in RMS and CVE for succesive temporal windows for each condition and subject

% Simulate Gaussian noise to create thresholds
Gaussian_noise = wgn(1000000,2,0);
window             = 1000; overlap = 500;  dts     = 100;  step     = window - overlap;
win_i = 1;
swin = (win_i*step/dts+1):(win_i*step/dts+window/dts+1);
dim                   = size(Gaussian_noise); tpoints = dim(1);
W                      = round((tpoints-window/dts-dts)/(window-overlap)*dts);
W;
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
noise_cve                 = zeros(dim(2),W);
for win_i = 0:W-1
    for ch_i = 1:dim(2)
        X = envelope(ch_i,(win_i*step/dts+1):(win_i*step/dts+window/dts+1));
        RMS_envelope(ch_i, win_i+1) = rms(X);
        noise_cve(ch_i, win_i+1) = std(X)/mean(X);
    end
end
mean(noise_cve(:));

% Compute shifts in RMS and CVE for succesive temporal windows for each condition and subject
Gaussian_CVE = sqrt((4-pi)/pi) % Gaussian CVE

switch threshold
    case 'Gaussian_Noise'
        thr_rhythm = Gaussian_CVE;
        thr_phasic = Gaussian_CVE;
    case 'Bandpass_Percentiles'
        thr_rhythm = prctile(noise_cve(:),10);
        thr_phasic = prctile(noise_cve(:),90);
end


for i_cond = 1 : length(list.conditions)
    for i_subj = 1:length(eval(['results.' list.conditions{i_cond }]));
        
        % RMS
        Y_a = RMS{i_cond, i_subj};
        for i = 1:length(Y_a)-1; dY_a(i) = Y_a(i) - Y_a(i+1); end
        Ex = dY_a < 0;
        %max_excursionLength = max(diff(find([1,diff(Ex),1]))); % empirical
        runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
        %Excursions = Ex(starts(find(runs == excursion_length)));
        Excursions = Ex(starts(runs > excursion_length));
        RMS_excursions{i_cond, i_subj} = length(Excursions);
        RMS_down{i_cond, i_subj}   = length(nonzeros(Excursions));
        RMS_up{i_cond, i_subj}   = length(nonzeros(~Excursions));
        
        %  CVE
        X_a = CVE{i_cond, i_subj};
        for i = 1:length(X_a)-1; dX_a(i) = X_a(i) - X_a(i+1); end
        Ex = dX_a < 0;
        %max_excursionLength = max(diff(find([1,diff(Ex),1])));
        runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
        %Excursions = Ex(starts(find(runs == excursion_length)));
        Excursions = Ex(starts(runs > excursion_length));
        CVE_excursions{i_cond, i_subj} = length(Excursions);
        CVE_down{i_cond, i_subj}  = length(nonzeros(Excursions));
        CVE_up{i_cond, i_subj}   = length(nonzeros(~Excursions));
        clear dX_a
        
        %  Departures towards rhythmic activity
        X_a = CVE{i_cond, i_subj};
        c = 0
        for i = 1:length(X_a)-1;
            if X_a(i) <  thr_rhythm
                c = c + 1
                dX_a(c) = X_a(i) - X_a(i+1);
            end
        end
        Ex = dX_a > 0;
        %max_excursionLength = max(diff(find([1,diff(Ex),1])));
        runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
        %Excursions = Ex(starts(find(runs == excursion_length)));
        Excursions = Ex(starts(runs > excursion_length));
        CVE_rhythm_departures{i_cond, i_subj} = length(Excursions);
        CVE_rhythm_up{i_cond, i_subj}  = length(nonzeros(Excursions));
        clear dX_a
        
        %  Departures towards phasic activity
        X_a = CVE{i_cond, i_subj};
        c = 0
        for i = 1:length(X_a)-1;
            if X_a(i) >  thr_phasic
                c = c + 1
                dX_a(c) = X_a(i) - X_a(i+1);
            end
        end
        Ex = dX_a > 0;
        %max_excursionLength = max(diff(find([1,diff(Ex),1])));
        runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
        %Excursions = Ex(starts(find(runs == excursion_length)));
        Excursions = Ex(starts(runs > excursion_length));
        CVE_phasic_departures{i_cond, i_subj} = length(Excursions);
        CVE_phasic_up{i_cond,i_subj}  = length(nonzeros(~Excursions));
        clear dX_a
    end
    RMS_dyn(i_cond,1,:) = (100/sum(vertcat(RMS_excursions{i_cond,:})))*sum(vertcat(RMS_down{i_cond,:}));
    RMS_dyn(i_cond,2,:) = (100/sum(vertcat(RMS_excursions{i_cond,:})))*sum(vertcat(RMS_up{i_cond,:}));
    CVE_dyn(i_cond,1,:) = (100/sum(vertcat(CVE_excursions{i_cond,:})))*sum(vertcat(CVE_down{i_cond,:}));
    CVE_dyn(i_cond,2,:) = (100/sum(vertcat(CVE_excursions{i_cond,:})))*sum(vertcat(CVE_up{i_cond,:}));
    CVE_departures(i_cond,1,:) = (100/sum(vertcat(CVE_rhythm_departures{i_cond,:})))*sum(vertcat(CVE_rhythm_up{i_cond,:}));
    CVE_departures(i_cond,2,:) = (100/sum(vertcat(CVE_phasic_departures{i_cond,:})))*sum(vertcat(CVE_phasic_up{i_cond,:}));
    
end

RMS_dyn
CVE_dyn
CVE_departures



