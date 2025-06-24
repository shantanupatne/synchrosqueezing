# Analysis of Synchrosqueezing Transforms for Non-Stationary Signals

This project explores and assesses SSTs as a powerful tool for the TF representation of non-stationary signals. SSTs leverage adaptive TF analysis by concentrating signal energy in localized regions of the TF plane, thereby providing an enhanced representation of time-varying signal characteristics. This project investigates the capabilities and limitations of SST in capturing the dynamic behavior of non-stationary signals, and the efficacy of its evolved versions, the adaptive SST and second-order SST. 

The goal of this project is to determine whether or not SSTs are a suitable TFR for non-stationary signals under a variety of conditions. The main conditions we will be considering are combined linear chirp signals, combined linear signals that cross over each other, and noisy variants of the two. This is to determine under which conditions an SST is most suitable for signal processing, and to determine when performing an adaptive or second-order transformation might be necessary when working with the SST.

---

1. Adaptive Array of Lienard-type Intermittent Oscillators (ALI) for instantaneous frequency estimation in high noise environment as described by [Time-frequency high-resolution for weak signal detection using chaotic intermittence (Costa et al., 2023)](https://doi.org/10.1016/j.dsp.2023.104160).

2. Adaptive STFT-based Synchrosqueezing Transform for noisy multicomponent signal TF representation as described by [Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation (Li et al., 2020)](https://doi.org/10.1016/j.sigpro.2019.07.024).

3. Adaptive CWT-based Synchrosqueezing Transform similar to above as described by [Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation (Li et al., 2020)](https://doi.org/10.1016/j.acha.2019.06.002).

Code for (1) obtained from [here](https://codeocean.com/capsule/9015585/tree/v1) 

MATLAB routines for (2) and (3) as available [here](http://www.cs.umsl.edu/~jiang/Jsoftware_SST.htm)

*Course Project for EEE505: Time-Frequency Signal Analysis.*
