# Rhythm Gaussianity

Here a minimal working example for revealing Rhythm Gaussianity by relating amplitude modulations to EEG signal energy.

* CVE = 0.523 reflects Gaussian noise
* CVE < 0.523 reflects rhythmic fluctuations (*e.g.* Kuramoto oscillations)
* CVE > 0.523 reflects phasic activity (*e.g.* avalanches)

> Hidalgo VM, Letelier JC, Diaz J, **The amplitude modulation pattern of Gaussian noise is a fingerprint of Gaussianity**, https://doi.org/10.48550/arXiv.2203.16253

## Introduction

Neuronal synchronization has long been assumed to be a main cause for the shaping of the alpha rhythm **[1]**. While some claimed it to be the footprint of synchronous coupled neuronal generators **[2]**, an idea akin to that of weakly coupled nonlinear oscillators **[3]**), others have noticed that the alpha rhythm can be related to a band-pass filter **[4]**. This idea was later corroborated by Motokawa, who, in a study that approximated the envelope of the alpha rhythm with a Rayleigh distribution **[5]**, demonstrated the alpha rhythm band-pass filter properties, and hence its relation to Gaussian noise **[6]**. Recently, a study by Hidalgo furthered this idea by proposing it as a method to assess the rhythm Gaussianity **[6]**. Here, we extend on this proposal by applying the central idea, that the coefficient of variation of the envelope can be used as a marker of rhythm Gaussianity, to study EEG alpha rhythm dynamics in two notoriously different experimental conditions. In what follows, we adapt this approach and extend its scope with different methodological constraints and experimental conditions. 

## **Example analysis**

## *Processing of EEG data*
Electrode signals from five occipital and parietal electrodes in N subjects in conditions A and B were re-referenced to the global average, normalized within channels, and band-pass filtered between 8-14 Hz using a zero-shift Butter-worth filter of order 4. The complex analytical signal was computed using the Hilbert transform applied to the whole time series length. Due to the removal of muscular artifacts, the coefficient of variation and the root median square (CVE and RMS) of the envelope had to be computed for sliding windows of 4 seconds of length and 25% overlap. Most artifacts occurred in condition A. The coefficient of variation of the envelope (CVE) was computed as the standard deviation divided by the mean of the envelope within each contiguous window (epoch). 

To visualize the RMS and the CVE distributions in each condition, the RMS envelope and the CVE values for each contiguous epoch were grouped across subjects and plotted as histograms (**Figure 1**). Subsequently, to visualize the relationship between the RMS and the CVE, we obtained bi-variate histograms as follows. First, the RMS values for each epoch were pooled across subjects and conditions to obtain minimum and maximum RMS values across conditions. Second, for each condition, the RMS and the CVE for each epoch were grouped into a two dimensional grid defined in CVE bins of size 0.005 (between 0 and 1) and RMS bins of 0.04 times the global maximum RMS envelope value (ranging from 0 to maximum). This step provided us with a common grid. Third, to visualize the resulting bi-variate distributions (the resulting bin counts for each each condition), we used contour plots and normalized them using the maximum bin count (**Figure 2**). 

<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/figures/Fig1.png" width=100%>

**Figure 1.** **(A)** Distribution of power in the alpha band (RMS). **(B)** Distribution of the coefficient of variation of the envelope (CVE). Light red corresponds to the condition A. Alpha power significantly decreased in the condition A (compared to that in the condition B; two-sample t-test, p-value <0.0001), whereas the CVE distributions were not significantly different between conditions (two-sample t-test, p-value = 0.0011). Data were grouped in a common space of 100 bins comprising the minimum and maximum values pooled over conditions.

  
### *Itinerancy in the RMS-CVE space* 
To quantify the itinerancy in the RMS-CVE space, the tendency for departures from Gaussianity, we proceeded as follows. First, for each subject and condition, we counted all sequences of contiguous epochs in which the RMS envelope and the CVE were increasing or decreasing. To do so, we defined a minimum sequence length (3 consecutive epochs, equivalent to 9 seconds) and identified the start of each sequence. Second, we summed the number of contiguous sequences with positive or negative shifts. Third, we computed, for each subject and condition, the percentage of the sequences.

To assess the departures from rhythm Gaussianity, we created a null distribution of CVE values by applying the procedure used to process the EEG data to Gaussian noise. This resulted in a distribution of CVE values with mean 0.523 (very close to the theoretical value for Gaussian noise: √(4 − π)/π ≈ 0.523). We then used the 10% and 90% percentile values of this null distribution as thresholds to define departures from Gaussianity. This allowed us to compute, for each subject and condition, the percentage of epochs in which the data departed from Gaussianity. 


<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/figures/Fig2.png" width=75%>

**Figure 2.** Contour plots illustrating the relationship between RMS envelope and CVE. Contours were obtained by linear interpolation of the bi-variate grid. Colorbar indicate bin count for each contour.
	
|   |   *RMS decreases*| *RMS increases*  | *Rhythmic CVE* | *Phasic CVE* | 
|---|---|---|---|---|
|  **A** |    50.7692 %|    49.2308 %|   1.1628 % | 3.4483 % |
|   **B**|    44.3011 %|    55.6989 %|   2.8302 % | 1.3333 % |

  
**Table 1.** Changes in the RMS and the CVE point to departures from rhythm Gaussianity. The first two columns indicate, for RMS, the percentage of contiguous epoch sequences longer than 9 seconds. The last two columns indicate, for CVE, the percentage of sequences below or above the threshold for Gaussian noise (the majority of epochs are between these values, hence the small percentages).


## Results
The mean of the RMS envelope distributions differed between the B and the A condition (two-sample t-test, p-value <0.0001). On the other hand, the mean of the CVE distributions did not (two-sample t-test, p-valjue = 0.0011). For sequences of contiguous epochs longer than 9 seconds, positive RMS shifts were increased in condition B as compared to condition A (**Table 1**, second column). In the case of CVE shifts, increased phasic activity was increased in condition B, as compared to conditon B (**Table 1**, last column). This indicates that departures from the Gaussianity in the direction of phasic activity were higher for condition A. On the other hand, departures towards oscillatory dynamics were increased for condition B (**Table 1**, third column). 
  
## Discussion and concluding remarks
Our analysis of itinerancy in the CVE space revealed, for conditon A, an increased tendency for departures from rhythm Gaussianity in the alpha band towards phasic activity. This result is of fundamental methodological importance, as it reveals, in a data driven manner, departures in the dynamic repertoire of rhythmic itinerancies, likely pointing to the disintegration of alpha oscillations (*e.g.* phase scattering, neuronal avalanches), hypothesized to occur more frequently in condition A **[7]**. Even though the observed increased departures from Gaussianity towards the direction of phasic activity (**Table 1**, last column) may indeed relate to changes in the mean EEG field caused by neuronal mechanisms **[8]**, it is likely that these too reflect the contribution of muscular artifacts. Finally, we pushed the limits of the method by defining epochs using a sliding window length of 4 seconds and 25%, instead of the 24 seconds prescribed in **[6]**. While this decision ensured that muscular artifacts could be efficiently rejected, we cannot discard the potentially limiting influence of this decision on the CVE estimation. In the future, these limitations could be addressed to improve the potential of this method to assess task-dependent dynamic repertoires in EEG rhythms.


## References

1. Singer, Wolf. (1999).**Neuronal synchrony: A versatile code for the definition of relations?** Neuron, vol. 24, no. 1, pp. 49–65.

2. Strogatz, Steven. H. (1994). **Norbert Wiener’s brain waves**. Frontiers in Mathematical Biology Berlin, Heidelberg: Springer, pp. 122–138.

3. Kuramoto, Yoshiki. (1975). **Self-entrainment of a population of coupled non-linear oscillators.**  International Symposium on Mathematical Problems in Theoretical Physics. Lecture Notes in Physics, vol 39. Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0013365 

4. Da Silva, F. H. L. et al. (1973). **Organization of thalamic and cortical alpha rhythms: Spectra and coherences**. Electroencephalography Clin. Neuriophysiol., vol. 35, no. 6, pp. 627–639.

5. Motokawa, Kōiti And Tosisada Mita. (1942). **Das Wahrscheinlichkeitsprinzip über Die Gehirn-Elektrische Erscheinungen Des Menschen**. Japanese Journal of Medical Sciences III (Biophysics 8): 63–67.

6. Hidalgo Victor M, Diaz J, Mpodozis J, Letelier JC. (2023). **Envelope Analysis of the Human Alpha Rhythm Reveals EEG Gaussianity**. IEEE Trans Biomed Eng. 2023 Apr;70(4):1242-1251. https://doi.org/10.1109/TBME.2022.3213840  

7. Jevri Hanna, Cora Kim, Stefan Rampp, Michael Buchfelder, Nadia Müller-Voggel. (2024) **Decreasing alpha flow releases task-specific processing paths**. Imaging Neuroscience 2024; 2 1–24. doi: https://doi.org/10.1162/imag_a_00117

8. Freeman, Walter J. (2000). **How Brains Make Up Their Minds**. Columbia University Press.

