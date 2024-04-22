# Rhythm Gaussianity

Here a minimal working example for revealing Rhythm Gaussianity by relating amplitude modulations to signal energy.

* CVE = 0.523 reflects Gaussian noise
* CVE < 0.523 reflects rhythmic fluctuations (e.g. Kuramoto oscillations)
* CVE > 0.523 reflects phasic activity (e.g. avalanches)



**References**

1. Hidalgo VM, Diaz J, Mpodozis J, Letelier JC. **Envelope Analysis of the Human Alpha Rhythm Reveals EEG Gaussianity**. IEEE Trans Biomed Eng. 2023 Apr;70(4):1242-1251. https://doi.org/10.1109/TBME.2022.3213840  

2. Hidalgo VM, Letelier JC, Diaz J. **The amplitude modulation pattern of Gaussian noise is a fingerprint of Gaussianity**. https://doi.org/10.48550/arXiv.2203.16253



## **Example analysis**

## *Preprocessing of EEG data*
Electrode signals from five occipital and parietal electrodes in N subjects in B and A conditions were re-referenced to the global average, normalized within channels, and band-pass filtered between 8-14 Hz using a zero-shift Butter-worth filter of order 4. The complex analytical signal was computed using the Hilbert transform applied to the whole time series length. Due to the removal of muscular artifacts, the coefficient of variation and the root median square (CVE and RMS) of the envelope had to be computed for sliding windows of 4 seconds of length and 25% overlap, instead of the 24 seconds prescribed in **[1]**. Most artifacts occurred in condition A. The coefficient of variation of the envelope (CVE) was computed as the standard deviation divided by the mean of the envelope within each contiguous window. 

To visualize the RMS and CVE distributions in each condition, the RMS envelope and the CVE values for each contiguous window were grouped across subjects and plotted as histograms (**Figure 1**). Subsequently, to visualize the relationship between RMS envelope and CVE, we obtained bivariate as follows. First, the RMS envelope values for each window were pooled across subjects and conditions to obtain minimun and maximum RMS values across conditions. Second, for each condition, the RMS envelope and CVE values for each window were grouped into a two dimensional grid defined in CVE bins of size 0.005 (between 0 and 1) and RMS envelope bins of 0.04 times the global maximum RMS envelope value (ranging from 0 to maximum). This step providesd us with a common grid. Third, to visualize the resulting bivariate distributions (the resulting bin counts for each each condition), we used contour plots and normalized them using the maximum bin count (**Figure 2**). 

<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/figures/Fig1.png" width=100%>

**Figure 1.** **(A)** Distribution of power in the alpha band (RMS envelope). **(B)** Distribution of the coefficient of variation of the envelope. Light red corresponds to the condition A. Alpha power significantly decreased in the condition A (compared to that in the condition B; two-sampled t-test, p-val <0.0001), whereas the CVE distribution were not significantly different bwteen conditions (two-sampled t-test, p-val = 0.0011). 

### *Itinerancies in the RMS-CVE space*  
To asses increases and decrease in RMS envelope and CVE throughout the recordings, we quantified sequences of contiguous windows in which the RMS envelope and CVE values were increasing or decreasing. To do so, we defined a minimum sequence length (3 consecutive windows, equivalent to 9 seconds) and identified the start of each sequence. Second, we summed the number of contiguous sequences with positive or negative shifts longer than 9 seconds. Third, we computed the percentage of sequences of contiguous shifts above the minimum sequence length. Finally, to quantify the dynamic CVE range, we proceeded as follows. First, we pooled CVE shifts from each condition to obtain a common pool. This gave us a symmetric distribution centered at zero. Second, we obtained the 5% and 95% percentile values of this common distribution and used them to define low and high CVE thresholds coomon for both cobnditions. Subtracting and adding these values to  the CVE of Gaussian noise (√(4 − π)/π ≈ 0.523) allowed us to compute, for each condition, the number of windows above and below these thresholds. Because CVE values below or above 0.523 are indicative rhythmic and phasic dynamics, this approach helped us define the dynamic range of the CVE in each condition.

<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/figures/Fig2.png" width=75%>

**Figure 2.** Contour plots illustrating the relationship between RMS envelope and CVE. Contours were obtained by linear interpolation of the bivariate grid. Colorbar indicate bin count for each contour.
	
            |   |   **A**| **B**  | 
            |---|---|---|
            |   RMS increases|    39.6825 %|    61.2613 %|  
            |   RMS decreases|    60.3175 %|     38.7387 %|   
            |   CVE increases|   52.1739 %|   48.3871 %|  
            |   CVE decreases|    47.8261 %|   51.6129 %|
  
**Table 1.** Changes in RMS envelope and CVE in point to departures from rhythm Gaussianity. 


### *Results*
The mean of the RMS envelope distributions differed between the B and the A condition (two-sampled t-test, p-val <0.0001). On the other hand, the mean of the CVE distributions did not (two-sampled t-test, p-val = 0.0011). For sequences of contiguous windows above 9 seconds, positive RMS envelope shifts were increased in the B condition. Events longer than 9 seconds in which alpha power decreased were 60.3175% for the A condition and 38.7387% for the B condition. In the other hand, alpha power increasing events were 39.6825% for the A condition and 61.2613% for condition B. These results indicate that alpha power increases was somewhat suppressed for condition A. In the case of CVE, results told a different story. For sequences of contiguous windows above 9 seconds, positive CVE shifts were increased in the A condition. Events longer than 9 seconds in which CVE decreased were 47.8261% for the A condition and 51.6129% for the B condition. In the other hand, CVE increasing events were 52.1739% for condition A and 48.3871% for condition B (**Table 1**). This indicates that departures from the alpha rhythm Gaussianity in the direction of phasic activity were higher for condition A, possibly pointing to the disintegration of oscilatory dynamics (e.g. phase scatering), increased neuronal avalanches or another related mechanism. On the other hand, departures towards oscillatory behavior were increased for condition B. Finally, state occupancy revealed a reduced repertoire for A, with occupancy values of 1.1788 % and 2.3772% of the total window count for rhythmic and phasic activity, respectively, compared to 2.1998% and 3.593% in the B condition, concomitant with increases in apha power. 
  
### *Discussion and concluding remarks*
Regardless of the reduced dynamic range in RMS envelope and CVE (revealed by our state occupancy analysis), our analysis of RMS envelope sequence shifts confirmed decreases in alpha power for the A condition. Moreover, for condition A, CVE sequence shifts revealed an increased tendency for departures from Gaussianity towards phasic activity. This last result is of fundamental methodological importance, as it reveals, in a data driven manner, departures in the dynamic repertoire of rhythmic itineraries, likely pointing to mechanisms such as neuronal avalanches, hypothesized to occur more frequently for condition A. Even though the observed increased departures from Gaussianity towards the direction of phasic activity in condition A (**Table 1**), CVE increases in condition A may relate to changes in the mean EEG field caused by neuronal mechanisms, it is likely that they too reflect the contribution of muscular artifacts. Moreover, the calculation of the sequences was done on time series concatenated across subjects, leaving open the possibility that windows belonging to different subjects appeared contiguous in our analysis. While this procedure may indeed introduce jumps and these jumps may contribute to apparent phasic dynamics, we believe that this effect is negligible, as the the minimum window length would in theory reduce this. Finally, we pushed the limits of the method by using a window length of 4 seconds. While this decision ensured that no artifacts were included in the data, we cannot discard the potentially limiting influence on the CVE estimation. In the future, the limitations here described could be assessed to improve this potential so this method can be further validated. To conclude, we acknowledge the potential of using rhythm Gaussianity analysis to assess dynamic repertoires in rhythmic EEG data variety of conditions in health and disease. 

