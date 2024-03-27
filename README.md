# rhythm-gaussianity

Revealing Rhythm Gaussianity by using Coefficient of Variation of the Envelope (CVE).

### Compute Coefficient of Variation (CVE) and Power (RMS)
* CVE = 0.523 reflects Gaussian noise
* CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
* CVE > 0.523 reflects phasic activity (e.g. avalanches)



* Reference:
 Hidalgo VM, Diaz J, Mpodozis J, Letelier JC. Envelope Analysis of the Human Alpha Rhythm Reveals EEG Gaussianity. IEEE Trans Biomed Eng. 2023 Apr;70(4):1242-1251. doi: 10.1109/TBME.2022.3213840. Epub 2023 Mar 21. PMID: 36223351.

### Here we provide a minimal working example in Matlab:

1. Load time series (as 2D array: *channels x time points*)
2. Normalize acorss channels (z-score)
3. Concatenate the time series in each channels (to a vector)
4. Filter all channels (as vector)
5. Obtain analytical signal: apply Hilbert transform
6. De-concatenate channels (back to 2D array)
7. Run sliding window to obtain, for each window:
    * Envelope: Compute the Coefficient of Variation (CVE) and the RMS envelope 
    * Phase: Compute local and global phase coherences using Kuramoto defintion (as *channel i to all* and each *across all channels*, respectively) 

Example result for a single subject without referencing:
  
<img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/rmsenvDyn_noref_freq_1.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/phaseDyn_noref_freq_1.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/ph-rmsenvDyn_Laplace_freq_1.png" width=30%>
<img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/rmsenvDyn_noref_freq_2.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/phaseDyn_noref_freq_2.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/ph-rmsenvDyn_Laplace_freq_2.png" width=30%>
<img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/rmsenvDyn_noref_freq_3.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/phaseDyn_noref_freq_3.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/ph-rmsenvDyn_Laplace_freq_3.png" width=30%>
  

    
* Example result for a single subject Laplace referencing:
  
<img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/rmsenvDyn_Laplace_freq_1.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/phaseDyn_Laplace_freq_1.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/ph-rmsenvDyn_Laplace_freq_1.png" width=30%>
<img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/rmsenvDyn_Laplace_freq_2.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/phaseDyn_Laplace_freq_2.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/ph-rmsenvDyn_Laplace_freq_2.png" width=30%>
<img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/rmsenvDyn_Laplace_freq_3.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/phaseDyn_Laplace_freq_3.png" width=30%><img src="https://github.com/nicogravel/rhythm-gaussianity/blob/main/mwe/ph-rmsenvDyn_Laplace_freq_3.png" width=30%>

