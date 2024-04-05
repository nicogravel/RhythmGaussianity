%% Compute Coefficient of Variation (CVE) and Power (RMS)
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
%
% 1. Load time series as 2D array: *channels x time points*
% 2. Normalize acorss channels using z-score
% 3. Concatenate the time series in each channel to a vector
% 4. Filter normalized time series
% 5. Obtain analytical signal: apply Hilbert transform
% 6. De-concatenate the analytical signal back to a 2D array with shape *channels x time points*
% 7. Run a sliding window to obtain, for each window:
%   * Envelopes: To then compute the Coefficient of Variation (CVE) and the RMS envelope. 
%   * Phases: To then compute local and global phase coherences as the Kuramoto order prameter for *channel i to all* and *across all channels*

clear all
close all

%% Band-pass filter
f_sampling = 250; % 1kHz
freqs_low  = [1, 8, 30]; %2 % Hz
freqs_high = [6, 14, 60]; %12; % Hz   (Hidalgo et al. use 45 Hz)
conditions = ['A', 'B'];

load('MWE.mat')

referencing = 2;
switch referencing
    case 1
        data      = data_unref
    case 2
        data      = data_ref
end

% Determine sliding window
tpoints = 49348;
seconds = tpoints/f_sampling
minutes = seconds/60
window   = 4000; overlap  = 2000;  dts   = 10;  % Windows length: 401 points
window   = 6000; overlap = 2000; dts     = 500;  % Windows length: 13 points
step     = window - overlap;
step     = window - overlap;
W        = round((tpoints-window/dts-dts)/(window-overlap)*dts);
W
win_i      = W; % examine windows, overlap and length!
win_length = (win_i*step/dts+1):(win_i*step/dts+window/dts+1);
size(win_length,2)

tic
% Loop over frequency 
for freq_i = 1:3
    
    
    %% Power Dynamics
    f_low  = freqs_low(freq_i);
    f_high = freqs_high(freq_i);
    
    f_sampling = 250; % 1kHz
    f_nrm_low   = f_low /(f_sampling/2);
    f_nrm_high  = f_high /(f_sampling/2);
    % Determine filter coefficients:
    [z,p,k] = butter(4,[f_nrm_low f_nrm_high],'bandpass');
    % Convert to zero-pole-gain filter parameter (recommended)
    sos = zp2sos(z,p,k);
    
    % Loop over conditions   
    for cond_i = 1:2
       
        
        %% Compute CVE
        
        % 1) Load time series (as 2D array)
        ts = data{cond_i}; size(ts);
        
        % 1.1) Check!
        %figure,plot(ts(:,100:200)')
        
        % 2)  Normalize across channels
        ts_norm = zscore(ts); size(ts_norm); % load and zscore data
        % 2.1) Check!
        %figure,plot(ts_norm(:,100:200)')
        
        
        % 3) Concatenate channels (to vector)
        concts = zeros(1,size(ts,1)*size(ts,2)); size(concts);
        for ch_i = 1:size(ts,1)
            id = size(ts,2)*ch_i;
            concts(1,(id-size(ts,2)+1):id) = ts_norm(ch_i,:);
        end
        size(concts);
        clear ts_norm
        % 3.1) Check!
        %figure,plot(concts(100:5000))
        
        % 4) Filter all channels (as vector)
        sig_flt = sosfilt(sos,concts); size(sig_flt); % apply filter
        clear concts
        
        % 4.1) Check!
        %figure,plot(sig_flt(100:1000))
        
        % 5) Obtain analytical signal: apply Hilbert transform
        hbert = hilbert(sig_flt); size(hbert);
        envel = abs(real(hbert));
        phase = angle(hbert);
        clear sig_flt hbert
        
        % 5.1) Check!
        %figure,plot(envel(100:1000))
        
        % 6) De-concatenate channels (back to 2D array)
        env_ch   = zeros(size(ts,1),size(ts,2)); size(env_ch);
        phi_ch   = zeros(size(ts,1),size(ts,2)); size(phi_ch);
        
        for ch_i = 1:size(ts,1)
            id = size(ts,2)*ch_i;
            env_ch(ch_i,:)   = envel((id-size(ts,2)+1):id);
            phi_ch(ch_i,:)   = phase((id-size(ts,2)+1):id);
        end
        size(env_ch);
        clear ts
        
        % 6.1) Check things are good!
        %figure,plot(env_ch(:,100:150)')
        %figure,plot(phi_ch(:,100:150)')
        
        
        % 7) Run sliding window
        tpoints = size(env_ch,2);
        step    = window - overlap;
        W       = round((tpoints-window/dts-dts)/(window-overlap)*dts);
        W
        for win_i = 0:W-1
            for ch_i = 1:size(env_ch,1)
                
                % Obtain windowed envelope
                wenv = env_ch(ch_i,(win_i*step/dts+1):(win_i*step/dts+window/dts+1));
                wenv_ch(ch_i, win_i+1, :) = wenv(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10)));
                
                % Obtain windowed phase
                wphi = phi_ch(ch_i,(win_i*step/dts+1):(win_i*step/dts+window/dts+1));
                wphi_ch(ch_i, win_i+1, :) = wphi(round(length(wenv)/10):(length(wenv)-round(length(wenv)/10)));
                
            end
        end
               
        % 8) Compute RMS envelope as function and CVE
        size(wenv_ch)
        env_all = reshape(wenv_ch,[size(wenv_ch,1)*size(wenv_ch,2),size(wenv_ch,3)]); size(env_all);
        cve  = std(env_all, 0,2)./mean(env_all,2); size(cve);
        pow  = rms(env_all,2); size(pow);
        
        % 9) Compute inter-link phase coherence for each channel and window
        for win_i = 0:W-1
            for ch_i = 1:size(wphi_ch,1)
                for tp_i = 1:size(wphi_ch,3)
                    dPhi_t(:,tp_i) = wphi_ch(ch_i, win_i+1, tp_i) - wphi_ch(:, win_i+1, tp_i);
                end
                % Compute interlink phase coherence (local Kuramoto order)
                for tp_i = 1:size(wphi_ch,3)
                    D = nonzeros(dPhi_t(:,tp_i)); % ignore self loops
                    D_t(tp_i) = abs(real(mean(exp(1i*wrapToPi(D)))));
                end
                mtst_ch(ch_i,win_i+1) = std(D_t);
            end
        end
        size(mtst_ch);
        mtst_all = reshape(mtst_ch,[size(mtst_ch,1)*size(mtst_ch,2)],1); size(mtst_all);
        clear dPhi_t
        
        % 10) Compute global phase coherence for each window
        for win_i = 0:W-1
            for ch_i = 1:size(wphi_ch,1)
                for ch_j = 1:size(wphi_ch,1)
                    for tp_i = 1:size(wphi_ch,3)
                        dPhi_t(ch_i,ch_j,tp_i) = wphi_ch(ch_i, win_i+1, tp_i) - wphi_ch(ch_j, win_i+1, tp_i);
                    end
                end
            end
            % Compute global phase coherence (global Kuramoto order)
            for tp_i = 1:size(wphi_ch,3)
                D = dPhi_t(:,:,tp_i);
                D = D(find(tril(ones(size(D)),-1)));
                D_t(tp_i) = abs(real(mean(exp(1i*wrapToPi(D))))); % ignore self loops
            end
            mtst_global(ch_i,win_i+1) = std(D_t);
        end
        size(mtst_global);
        clear dPhi_t
        
        % Save results
        results.CVE{freq_i, cond_i}         = cve;         % This you plot against RM envelope
        results.POW{freq_i, cond_i}         = pow;
        results.MTST_ch{freq_i, cond_i}     = mtst_all;    % This you can plot against RMS envelope or CVE
        results.MTST_global{freq_i, cond_i} = mtst_global; % This you can average across all windows
        
        clear cve pow env_all mtst_global mtst_ch mtst_all wenv_ch wphi_ch
             
    end   
end
toc

save test_results results

