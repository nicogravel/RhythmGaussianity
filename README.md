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
Electrode signals from five occipital and parietal electrodes in 20 subjects in B and A conditions were re-referenced to the global average, normalized within channels, and band-pass filtered between 8-14 Hz using a zero-shift Butter-worth filter of order 4. The complex analytical signal was computed using the Hilbert transform applied to the whole time series length. The coefficient of variation and the root median square (CVE and RMS) of the envelope were computed for sliding windows of 4 seconds length and 25% overlap. No artifacts rejection was applied. The coefficient of variation of the envelope (CVE) was computed as the standard deviation divided by the mean of the envelope within each contiguous window. 

### *RMS envelope as a function of CVE*  
For each condition, the RMS envelope and the CVE  values for each contiguous window were grouped across subjects, binned (N=100) and plotted as histograms (Figure 1). To compare the RMS envelope and CVE distributions between conditions, we used two-sampled t-test. Subsequently, to visualize the relationship between RMS envelope and the CVE in each condition, bi-variate histograms were obtained as follows. First, the RMS and CVE values for each window were pooled across subjects and the two conditions to obtain global minimum and the maximum values for each estimate. Second, for each condition, the RMS envelope and CVE from each window were grouped into the cells of a two-dimensional grid defined in CVE bins of 0.005 between 0 and 1 and RMS envelope bins of 0.04 times the global maximum CVE between 0 and the global maximum value of the CVE (obtained from both conditions as explained previously). This step provides with a common histogram grid for both condition. Third, to compare the resulting bi-variate distributions (the resulting bin counts in each each condition), we used contour plots and normalized them using the maximum bin count (Figure 2). 

<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/figures/Fig1.png" width=100%>

**Figure 1.** (A) Distribution of power in the alpha band (RMS envelope). The EEG signal for  was  (B) Distribution of the coefficient of variation of the envelope. Alpha power significantly decreased in the A, compared to that of the B condition, whereas the CVE distribution were not significantly different. Light red corresponds to the A condition. 

### *Shifts in RMS and CVE*  
To asses the tendency in each condition to decrease or increase power (RMS) and CVE throughout the recordings, we quantified sequences of shifts in the RMS sign. First, we defined a minimum sequence length (3 consecutive windows, equivalent to 9 seconds) and identified the start of each sequence. Second, we summed the number of sequences with positive or negative shifts. Third, we computed the percentage of signed sequence shifts above the minimum sequence length.  
 
### *State occupancy* 
To quantify the range of dynamic activity in each conditions captured by the CVE, we proceeded as follows. First, we computed the global CVE shifts (all window shifts pooled for each condition). This gave us a symmetric distribution centered at zero.  Second, we obtained the 5% and 95% percentile values of this distribution and used them to define low and high CVE thresholds by subtracting or adding them (respectively) to the CVE of Gaussian noise (√(4 − π)/π ≈ 0.523). This allowed us to compute, for each condition, the number of temporal windows occupying the lower or higher ends of the global CVE space. Values below 0.523 in CVE indicate rhythmic activity whereas higher values indicate phasic activity.

<img src="https://github.com/nicogravel/RhythmGaussianity/blob/main/figures/Fig2.png" width=75%>

**Figure 2.** Contour plots reflecting the relationship between EEG amplitude and CVE.  Distribution of power in the alpha band (RMS envelope in the y-axis) as a function of the coefficient of variation of the envelope (CVE in the x-axis) for temporal windows pooled across subjects. Histogram bin counts were made comparable between the A (A) and B condition (B) by normalizing the bi-variate grid and setting the colormap to the maximum bin count. Contours were obtained by linear interpolation of the bi-variate grid. Colorbar indicate bin count for each contour.
	
|   |   **A**| **B**  | 
|---|---|---|
|   RMS increases|    39.6825 %|    61.2613 %|  
|   RMS decreases|    60.3175 %|     38.7387 %|   
|   CVE increases|   52.1739 %|   48.3871 %|  
|   CVE decreases|    47.8261 %|   51.6129 %|
**Table 1.** Changes in RMS envelope and CVE in point to departures from rhythm Gaussianity. To quantify the percentage of increases in RMS envelope and CVE departures n


### *Results and concluding remarks*
The mean of the RMS envelope distributions differed between the B and the A condition (two-sampled t-test, p-val <0.0001). On the other hand, the mean of the CVE distributions did not (two-sampled t-test, p-val = 0.0011). For sequences of contiguous windows above 9 seconds, positive RMS shifts were increased in the B condition. Events longer than 9 seconds in which alpha power decreased were 60.3175% for the A condition and 38.7387% for the B condition. In the other hand, alpha power increasing events were 39.6825% for the A condition and 61.2613% for the B condition. These results indicate that alpha power increases was somewhat suppressed for A.  In the case of CVE, results told a different story. For sequences of contiguous windows above 9 seconds, positive CVE shifts were increased in the A condition. Events longer than 9 seconds in which CVE decreased were 47.8261% for the A condition and 51.6129% for the B condition. In the other hand, CVE increasing events were 52.1739% for the A condition and 48.3871% for the B condition. This indicates that departures from the alpha rhythm Gaussianity in the direction of phasic activity were higher in the A condition, possibly pointing to increased neuronal avalanches. On the other hand, departures towards oscillatory behavior were increased in the B conditions. Finally, state occupancy revealed a reduced repertoire for A, with occupancy values of 1.1788 % and 2.3772% of the total window count for rhythmic and phasic activity, respectively, compared to 2.1998% and 3.593% in the B condition. Regardless of the reduced dynamic range in RMS and CVE (revealed by our state occupancy analysis), our analysis of RMS sequence shifts confirmed decreases in alpha power for the A condition. Moreover, in the A condition, CVE sequence shifts revealed increased tendency of event departing from Gaussianity towards phasic activity. This last result is of fundamental methodological importance, as it reveals, in a data driven manner, departures in the dynamic repertoire of rhythmic itineraries, likely pointing to mechanisms such as neuronal avalanches, hypothesized to occur more frequently in the psychedelic state. 

