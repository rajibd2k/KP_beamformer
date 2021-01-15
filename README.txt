R. Sharma, I. Cohen and J. Benesty, “Adaptive and hybrid Kronecker product beamforming for far-field speech signals,” Speech Communication, Elsevier, Vol. 120, pp. 42-52, 2020. 
https://doi.org/10.1016/j.specom.2020.04.001.



For synthetic speech
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1. Run Database_Synthetic_a.m

2. Run Database_Frames.m

3. Run Statistics_script_0.m

4. Run Filters_script_0.m, Filters_script_1.m

5. Run Metrics_iterations_M2.m, Metrics_snapshots.m, Metrics_iSIR.m, Metrics_PowerPattern.m




For real speech
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
One needs to have the TIMIT database in the system, resampled to 8 kHz.

1. Run Database.m 

2. Run Data.m

3. Filters_script_b.m

4. Run the scripts for generating the metrics and the figures.
