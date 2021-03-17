Most of these notebooks are only interesting if you have the original data. They are meant to be able to repeat the experiments and to calculate large amounts of directionalities at once, instead of having to go through the pipeline step-by-step for each contact pair


```run_all_computational_time```


```avg_frag_rmse_check```
Contains the code that is used to calculate the RMSE of average fragment.
* Experiment 1: The average RMSE of the kabsch alignment is calculated.
* Experiment 2: The RMSE if the labels are averaged.
* Experiment 3: The end RMSE, where k-means was used to reset the labels if neccessary (RMSE > 0.1).


```calc_all_directionalities```
Contains the code that is used to make directionality tables for all contact pairs. A result would look like the following table.
![Directionality table](../../results/plots/directionalities_10_03_rcome_ret_kmeans_res05_free_volume.svg)

Density calculations need to be already performed, and the volumes of the contact pairs need also already to be calculated.


```compare_rmses_alignment```
In this notebook, we find out what the best way is to align the fragments. We do this by comparing different rmse's of different superimposition (alignment) methods.

One of the results is the following:
![RMSEs of different alignments](../../results/plots/comparing_rmse_kabsch_rotation_rc6h5_r2co.svg)


```compression_test```
TODO
This notebook not only plots already existing results, but also runs the compression algorithm itself. 

bla bla

compression test example


```volume_central_groups```
![Volume errors](../../results/plots/volumes_error.svg)


```volume_R_groups```
Reasoning about why R groups should be half in the min volume, and not in the max volume. This is done with pictures like the following:
![R groups](../../results/plots/overlap_R_volumes_RNO2.png)

RNO2
with R (total)                : 248.6
with half R                   : 140.3      43.55 % difference
with R in min, not in max     : 153.7      38.17 % difference
with half R in min, not in max: 163.6      34.19 % difference

RCOMe
with R (total)                : 264.8
with half R                   : 156.5      40.89 % difference
with R in min, not in max     : 181.8      31.33 % difference
with half R in min, not in max: 190.1      28.21 % difference

RC6H5
with R (total)                : 387.9
with half R                   : 233.6      39.78 % difference
with R in min, not in max     : 247.5      36.21 % difference
with half R in min, not in max: 258.2      33.45 % difference

RC6F5
with R (total)                : 498.5
with half R                   : 305.0      38.82 % difference
with R in min, not in max     : 322.3      35.36 % difference
with half R in min, not in max: 335.2      32.76 % difference


```directionality_dependency_res_and_threshold```



```plot_fragments_with_labels```


