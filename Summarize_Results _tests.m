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


list.conditions = {'Mushroom', 'Placebo'};

%% Load results
for condition = 1 : length(list.conditions)
    fprintf('\n Printing %s:\n', list.conditions{condition})
    data = eval(['results.' list.conditions{condition}]);
    for file = 1: length(data)
        for i_chan = 1:5
            CVE_{condition,file,i_chan}  =  data{file}.cve(i_chan,:);  % CVE
            RMS_{condition,file,i_chan}  =  data{file}.pow(i_chan,:); % RMS envelope
        end
    end
end

%% Compute differences in RMS and CVE
RMS{1} = squeeze(vertcat(RMS_(1,1:size(results.Mushroom,2),:)));
RMS{2} = squeeze(vertcat(RMS_(2,1:size(results.Placebo,2),:)));
CVE{1} = squeeze(vertcat(CVE_(1,1:size(results.Mushroom,2),:)))
CVE{2} = squeeze(vertcat(CVE_(2,1:size(results.Placebo,2),:)));

% RMS
[H0,pval] = kstest2(horzcat(RMS{1}{:,:}),horzcat(RMS{2}{:,:}))  % CVE  differences are not significative

% CVE
[H0,pval] = ttest2(horzcat(CVE{1}{:,:}),horzcat(CVE{2}{:,:}),'tail','left','alpha',0.001)  % CVE  differences are not significative


%% CVE histograms
figure,
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-600 pos(2) 600, 250]); % Set plot size
set(gcf, 'color', 'w'); % Set figure background
Xa =  horzcat(RMS{1}{:,:})';
Xb =  horzcat(RMS{2}{:,:})';
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
Xa =  horzcat(CVE{1}{:,:})';
Xb =  horzcat(CVE{2}{:,:})';
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
X_1 =  horzcat(CVE{1}{:,:})'; %Y_a(:); %vertcat(CVE{1,:}); % CVE
Y_1 =  horzcat(RMS{1}{:,:})'; % X_a(:); %vertcat(RMS{1,:}); % RMS envelope
X_2 =  horzcat(CVE{2}{:,:})'; %Y_b(:); %vertcat(CVE{2,:}); % CVE
Y_2 =   horzcat(RMS{2}{:,:})'; %X_b(:); %vertcat(RMS{2,:});  % RMS envelope
minX  = min([X_1(:); X_2(:)]); maxX  = max([X_1(:); X_2(:)]);
minY  = min([Y_1(:); Y_2(:)]); maxY  = max([Y_1(:); Y_2(:)]);
binfactor = 25;
maxCount = 250;
for cond_i = 1:2
    subplot(1,2,cond_i);
    X =  horzcat(CVE{cond_i}{:,:})';
    Y =  horzcat(RMS{cond_i}{:,:})';
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
    ylim([minY 1.8]);
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
    clear dY dX dX_a dX_b
    for i_subj = 1:length(eval(['results.' list.conditions{i_cond }]))
        for i_chan = 1 : 5
            
            % RMS
            Y =RMS{i_cond}{i_subj,i_chan}; %RMS{i_cond, i_subj};
            for i = 1:length(Y)-1; dY(i) = Y(i) - Y(i+1); end
            Ex = dY < 0;
            %max_excursionLength = max(diff(find([1,diff(Ex),1]))); % empirical
            runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
            %Excursions = Ex(starts(find(runs == excursion_length)));
            Excursions = Ex(starts(runs > excursion_length));
            RMS_excursions{i_cond, i_subj,i_chan} = length(Excursions);
            RMS_down{i_cond, i_subj,i_chan}   = length(nonzeros(Excursions));
            RMS_up{i_cond, i_subj,i_chan}   = length(nonzeros(~Excursions));
            
            %  CVE
            X = CVE{i_cond}{i_subj,i_chan}; %CVE{i_cond, i_subj};
            for i = 1:length(X)-1; dX(i) = X(i) - X(i+1); end
            Ex = dX < 0;
            %max_excursionLength = max(diff(find([1,diff(Ex),1])));
            runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
            %Excursions = Ex(starts(find(runs == excursion_length)));
            Excursions = Ex(starts(runs > excursion_length));
            CVE_excursions{i_cond, i_subj,i_chan} = length(Excursions);
            CVE_down{i_cond, i_subj,i_chan}  = length(nonzeros(Excursions));
            CVE_up{i_cond, i_subj,i_chan}   = length(nonzeros(~Excursions));
            
            %  Departures towards rhythmic activity
            X = CVE{i_cond}{i_subj,i_chan}; %CVE{i_cond, i_subj};
            c = 0;
            for i = 1:length(X)-1;
                if X(i) <  thr_rhythm
                    c = c + 1;
                    dX_a(c) = X(i) - X(i+1);
                end
            end
            size(dX_a)
            Ex = dX_a > 0;
            %max_excursionLength = max(diff(find([1,diff(Ex),1])));
            runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
            %Excursions = Ex(starts(find(runs == excursion_length)));
            Excursions = Ex(starts(runs > excursion_length));
            CVE_rhythm_departures{i_cond, i_subj,i_chan} = length(Excursions);
            CVE_rhythm_up{i_cond, i_subj,i_chan}  = length(nonzeros(Excursions));
            
            %  Departures towards phasic activity
            X = CVE{i_cond}{i_subj,i_chan};  %CVE{i_cond, i_subj};
            c = 0;
            for i = 1:length(X)-1;
                if X(i) >  thr_phasic
                    c = c + 1;
                    dX_b(c) = X(i) - X(i+1);
                end
            end
            size(dX_b)
            Ex = dX_a > 0;
            %max_excursionLength = max(diff(find([1,diff(Ex),1])));
            runs = diff(find([1,diff(Ex),1])); starts = find([1,diff(Ex)]);
            %Excursions = Ex(starts(find(runs == excursion_length)));
            Excursions = Ex(starts(runs > excursion_length));
            CVE_phasic_departures{i_cond, i_subj, i_chan} = length(Excursions);
            CVE_phasic_up{i_cond,i_subj,i_chan}  = length(nonzeros(~Excursions));
        end
        
    end
end

%% Across subjects and channels
for i_cond = 1 : length(list.conditions)
    clear dY dX dX_a dX_b
    for i_subj = 1:length(eval(['results.' list.conditions{i_cond }]))
        RMS_dyn(i_cond,1,:) = (100/sum(vertcat(RMS_excursions{i_cond,:,:})))*sum(vertcat(RMS_down{i_cond,:,:}));
        RMS_dyn(i_cond,2,:) = (100/sum(vertcat(RMS_excursions{i_cond,:,:})))*sum(vertcat(RMS_up{i_cond,:,:}));
        CVE_dyn(i_cond,1,:) = (100/sum(vertcat(CVE_excursions{i_cond,:,:})))*sum(vertcat(CVE_down{i_cond,:,:}));
        CVE_dyn(i_cond,2,:) = (100/sum(vertcat(CVE_excursions{i_cond,:,:})))*sum(vertcat(CVE_up{i_cond,:,:}));
        CVE_departures(i_cond,1,:) = (100/sum(vertcat(CVE_rhythm_departures{i_cond,:,:})))*sum(vertcat(CVE_rhythm_up{i_cond,:,:}));
        CVE_departures(i_cond,2,:) = (100/sum(vertcat(CVE_phasic_departures{i_cond,:,:})))*sum(vertcat(CVE_phasic_up{i_cond,:,:}));
    end
end


RMS_dyn
CVE_dyn
CVE_departures


%% For individual channels
clear RMS_dyn CVE_dyn CVE_departures
for i_cond = 1 : length(list.conditions)
    clear dY dX dX_a dX_b
    for i_subj = 1:length(eval(['results.' list.conditions{i_cond }]))
        for i_chan = 1:5            
            RMS_dyn(i_cond,i_chan,1,:) = (100/sum(vertcat(RMS_excursions{i_cond,:,i_chan})))*sum(vertcat(RMS_down{i_cond,:,i_chan}));
            RMS_dyn(i_cond,i_chan,2,:) = (100/sum(vertcat(RMS_excursions{i_cond,:,i_chan})))*sum(vertcat(RMS_up{i_cond,:,i_chan}));
            CVE_dyn(i_cond,i_chan,1,:) = (100/sum(vertcat(CVE_excursions{i_cond,:,i_chan})))*sum(vertcat(CVE_down{i_cond,:,i_chan}));
            CVE_dyn(i_cond,i_chan,2,:) = (100/sum(vertcat(CVE_excursions{i_cond,:,i_chan})))*sum(vertcat(CVE_up{i_cond,:,i_chan}));
            CVE_departures(i_cond,i_chan,1,:) = (100/sum(vertcat(CVE_rhythm_departures{i_cond,:,i_chan})))*sum(vertcat(CVE_rhythm_up{i_cond,:,i_chan}));
            CVE_departures(i_cond,i_chan,2,:) = (100/sum(vertcat(CVE_phasic_departures{i_cond,:,i_chan})))*sum(vertcat(CVE_phasic_up{i_cond,:,i_chan}));
        end
        
    end
end

RMS_dyn
CVE_dyn
CVE_departures



%% Differences
for i_cond = 1 : length(list.conditions)
    clear dY dX dX_a dX_b
    for i_subj = 1:length(eval(['results.' list.conditions{i_cond }]))
        RMS_dyn(i_cond,1,:) = (100/sum(vertcat(RMS_excursions{i_cond,:,:})))*sum(vertcat(RMS_down{i_cond,:,:}));
        RMS_dyn(i_cond,2,:) = (100/sum(vertcat(RMS_excursions{i_cond,:,:})))*sum(vertcat(RMS_up{i_cond,:,:}));
        CVE_dyn(i_cond,1,:) = (100/sum(vertcat(CVE_excursions{i_cond,:,:})))*sum(vertcat(CVE_down{i_cond,:,:}));
        CVE_dyn(i_cond,2,:) = (100/sum(vertcat(CVE_excursions{i_cond,:,:})))*sum(vertcat(CVE_up{i_cond,:,:}));
        CVE_departures(i_cond,1,:) = (100/sum(vertcat(CVE_rhythm_departures{i_cond,:,:})))*sum(vertcat(CVE_rhythm_up{i_cond,:,:}));
        CVE_departures(i_cond,2,:) = (100/sum(vertcat(CVE_phasic_departures{i_cond,:,:})))*sum(vertcat(CVE_phasic_up{i_cond,:,:}));
    end
end


RMS_dyn
CVE_dyn
CVE_departures


%% Differences between conditions

% Differences in power attenuation
[H0,pval] = kstest2(vertcat(RMS_down{1,:,:}),vertcat(RMS_down{2,:,:}),'tail','smaller')  % Different

[H0,pval] = kstest2(vertcat(RMS_down{1,:,1:2}),vertcat(RMS_down{2,:,1:2}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(RMS_down{1,:,3}),vertcat(RMS_down{2,:,3}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(RMS_down{1,:,4:5}),vertcat(RMS_down{2,:,4:5}),'tail','smaller')  % Different



% Differences in power increase
[H0,pval] = kstest2(vertcat(RMS_up{1,:,:}),vertcat(RMS_up{2,:,:}),'tail','larger')           % Not different

[H0,pval] = kstest2(vertcat(RMS_up{1,:,1:2}),vertcat(RMS_up{2,:,1:2}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(RMS_up{1,:,3}),vertcat(RMS_up{2,:,3}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(RMS_up{1,:,4:5}),vertcat(RMS_up{2,:,4:5}),'tail','smaller')  % Not different


% Difference betwen conditions for rythm
[H0,pval] = kstest2(vertcat(CVE_down{1,:,:}),vertcat(CVE_down{2,:,:}),'tail','smaller')   % Not different


[H0,pval] = kstest2(vertcat(CVE_down{1,:,1:2}),vertcat(CVE_down{2,:,1:2}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(CVE_down{1,:,3}),vertcat(CVE_down{2,:,3}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(CVE_down{1,:,4:5}),vertcat(CVE_down{2,:,4:5}),'tail','smaller')  % Not different


% Difference betwen conditions for phase
[H0,pval] = kstest2(vertcat(CVE_up{1,:,:}),vertcat(CVE_up{2,:,:}),'tail','smaller')            % Different

[H0,pval] = kstest2(vertcat(CVE_up{1,:,1:2}),vertcat(CVE_up{2,:,1:2}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(CVE_up{1,:,3}),vertcat(CVE_up{2,:,3}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(CVE_up{1,:,4:5}),vertcat(CVE_up{2,:,4:5}),'tail','smaller')  % Not different


%% Differences within conditions

[H0,pval] = kstest2(vertcat(RMS_down{1,:,1:2}),vertcat(RMS_down{1,:,4:5}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(RMS_down{2,:,1:2}),vertcat(RMS_down{2,:,4:5}),'tail','smaller')  % Not different


[H0,pval] = kstest2(vertcat(RMS_up{1,:,1:2}),vertcat(RMS_up{1,:,4:5}),'tail','smaller')           % Not different
[H0,pval] = kstest2(vertcat(RMS_up{2,:,1:2}),vertcat(RMS_up{2,:,4:5}),'tail','smaller')           % Not different


[H0,pval] = kstest2(vertcat(CVE_down{1,:,1:2}),vertcat(CVE_down{1,:,4:5}),'tail','smaller')  % Not different
[H0,pval] = kstest2(vertcat(CVE_down{2,:,1:2}),vertcat(CVE_down{2,:,4:5}),'tail','smaller')  % Not different


[H0,pval] = kstest2(vertcat(CVE_up{1,:,1:2}),vertcat(CVE_up{1,:,4:5}),'tail','smaller')           % Not different
[H0,pval] = kstest2(vertcat(CVE_up{2,:,1:2}),vertcat(CVE_up{2,:,4:5}),'tail','smaller')           % Not different



%% Differences within conditions

[H0,pval] = kstest2(vertcat(RMS_down{1,:,:}),vertcat(RMS_down{1,:,:}))    % Not different

[H0,pval] = kstest2(vertcat(CVE_down{1,:,:}),vertcat(CVE_up{1,:,:}))     % Not different


[H0,pval] = kstest2(vertcat(RMS_down{2,:,:}),vertcat(RMS_up{2,:,:}))   % Not different

[H0,pval] = kstest2(vertcat(CVE_down{2,:,:}),vertcat(CVE_up{2,:,:}))    % Not different
