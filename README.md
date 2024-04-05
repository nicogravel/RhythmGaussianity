# Rhythm Gaussianity

Revealing Rhythm Gaussianity by using the Coefficient of Variation of the Envelope (CVE).. **work in progres**...

* CVE = 0.523 reflects Gaussian noise
* CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
* CVE > 0.523 reflects phasic activity (e.g. avalanches)



**References**

Hidalgo VM, Diaz J, Mpodozis J, Letelier JC. Envelope Analysis of the Human Alpha Rhythm Reveals EEG Gaussianity. IEEE Trans Biomed Eng. 2023 Apr;70(4):1242-1251. doi: 10.1109/TBME.2022.3213840. Epub 2023 Mar 21. PMID: 36223351.  

Hidalgo VM, Letelier JC, Diaz J. The amplitude modulation pattern of Gaussian noise is a fingerprint of Gaussianity. https://doi.org/10.48550/arXiv.2203.16253

Here we provide two minimal working examples. A *toolbox free* example using Matlab (mwe folder, *stil exploratory*), and an working example using Matlab and Fieldtrip (mwe_fieldtrip folder).

The following plots were produced for EEG data for a single subject (sampling rate: 1000Hz). Electrode signals from occipital-parietal electrode were re-ferenceed to the average and band-pass filtered between 8-13 Hz (using a zero-shift Butterworth oforder 4). The complex analytical signal was computed using the Hilbert transform on the concatenated data. Subsequently, the CVE, RMS envelope and Kuramoto order parameters were computed for sliding windows (24 seconds long with 90% overlap. This must change to 50% when scaling up the analysis to compute and groupe results across subjects). No artefacts rejection was applied. In the following plots, the capitol leters A and B denote conditions. In parenthesis, the number of sliding windows for each condition.
 

**Normalisation across electrodes**   
<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/rmsenv-cve_alpha_across.png" width=50%><img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/kurvar-rms_alpha_across.png" width=50%>

**Normalisation within electrodes**  
<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/rmsenv-cve_alpha_within.png" width=50%><img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/kurvar-rms_alpha_within.png" width=50%>

**Distributions for CVE and the standard deviation of phase coherence**  
Gray and light red histograms corresponds to condition A and B respetively.  
<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/hist_cve_alpha.png" width=50%><img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/mwe_fieldtrip/kurvar_alpha.png" width=50%>
