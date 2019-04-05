# ICASSP_2019_Localized_random_sampling
Codes and notes of our ICASSP 2019 paper on CFO robust CS

# ICASSP_2019_Localized_random_sampling

IMPLEMENTATION DETAILS 

The file main.m in main_program folder contains the main program corresponding to [1].

Experiment 1 in [1]: 
Set resultid=1 (Line 26) in main.m and run.

In this experiment, the channel measurements obtained with an IID random phase shift-based CS matrix and the proposed design are perturbed by the same CFO and the same noise. Then, the perturbed measurements are fed to a CS algorithm (OMP) that outputs a sparse channel estimate. The channel estimate obtained with OMP may not be the same as the ideal one, as CFO induces phase errors in the samples. From a beam alignment perspective, a straightforward way to check the "closeness" of the estimate to the true channel is to compare the beam misalignments that arise with beamforming using the estimated solution. It can be noticed that CS with the proposed solution results in small beam misalignments when compared to CS with the commonly used random phase shift-based design.

Experiment 2 in [1]:
Set resultid=2 (Line 26) in main.m and run.

In this experiment, the SNR achieved with five different strategies is evaluated. They correspond to perfect channel information, CS with the proposed training + beam broadening (to account for the misalignments), CS with the common random training, CS with the common random training + beam broadening, and a non-coherent algorithm called Agile-Link.
The result in experiment 2 indicates that the proposed CS matrix design with beam broadening, is robust to large CFO errors over the one that uses the common random phase shift-based design. 

For more details about our research visit :- 

http://www.profheath.org/

https://sites.google.com/site/nitinmyers/home

[1] N. J. Myers, and R. W. Heath Jr., "Localized random sampling for robust compressive beam alignment", to appear in proc. of IEEE ICASSP 2019
