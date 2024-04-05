# Rhythm Gaussianity

Revealing Rhythm Gaussianity by using the Coefficient of Variation of the Envelope (CVE).. **work in progres**...

* CVE = 0.523 reflects Gaussian noise
* CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
* CVE > 0.523 reflects phasic activity (e.g. avalanches)



**References**

1. Hidalgo VM, Diaz J, Mpodozis J, Letelier JC. **Envelope Analysis of the Human Alpha Rhythm Reveals EEG Gaussianity**. IEEE Trans Biomed Eng. 2023 Apr;70(4):1242-1251. https://doi.org/10.1109/TBME.2022.3213840  

2. Hidalgo VM, Letelier JC, Diaz J. **The amplitude modulation pattern of Gaussian noise is a fingerprint of Gaussianity**. https://doi.org/10.48550/arXiv.2203.16253

Here we provide two minimal working examples. A *toolbox free* example using Matlab for illustrative purpose (mwe folder), and a working example using Matlab and the toolbox Fieldtrip (mwe_fieldtrip folder).

## **Example results**

The following plots were obtained from EEG data for a single subject (sampling rate: 1000Hz). Electrode signals from five occipital and parietal electrodes were re-referenced to the global average and band-pass filtered between 8-13 Hz using a zero-shift Butterworth oforder four. The complex analytical signal was computed using the Hilbert transform applied to the whole time series length. The CVE, RMS envelope and Kuramoto order parameters were computed for sliding windows of 24 seconds length and 90% overlap (this must change to 50% when scaling up the analysis). No artefacts rejection was applied. In the following plots, the capital leters A and B indicate arbitrary conditions. In parenthesis, the number of sliding windows for each condition (#).
 
### **RMS envelope as a function of CVE and phase coherence**  
  
For data normalised *across* electrodes:     
<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/rmsenv-cve_alpha_across.png" width=50%><img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/kurvar-rms_alpha_across.png" width=50%>

For data normalised *within* electrodes:    
<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/rmsenv-cve_alpha_within.png" width=50%><img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/kurvar-rms_alpha_within.png" width=50%>

### **Distributions for CVE and the standard deviation of phase coherence**  
  
Gray and light-red histograms corresponds to condition A and B:  
<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/hist_cve_alpha.png" width=50%><img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/kurvar_alpha.png" width=50%>
