**[V1 01/04/2024]** Here is the python code used in [ArXiv:2403.14580](https://arxiv.org/abs/2403.14580) to estimate the dipole and boost the mocks of the QSO eBOSS, LRG eBOSS, CMASS eBOSS, CMASS BOSS and LOWZ BOSS catalogs. In this version we are using quantil binning as an example, as it is faster, but you easily can change it to any binning strategy by modifing the "quantil_binning_number_list_s" function (used in the "z_dipole_estimator" function).


---

To **estimate the dipole** you should run the function: 
```
z_dipole_estimator(z_array,z_hat,w_syst,n_bins,n_threads=your_number_of_threads) 
```
The NGC and SGC (or any other sources with different monopoles) should be arranged in different rows in the z_array, z_hat and w_syst arrays.  

> z_array $\rightarrow$ is the observed 1+z array for each object
>
> z_hat $\rightarrow$ the array with the x, y and z components of the direction vector for each object
>
> w_syst $\rightarrow$ the systematic weight for each objetc.

You can change the grid (in beta) for all 3 iterations during the minimization of the estimator (default: x_min=-0.0021,x_max=0.0021,y_min=-0.0021,y_max=0.0021,z_min=-0.0021,z_max=0.0021,iter1_step=0.00007,iter2_step=0.000035,iter3_step=0.00001133333).


---

To **apply the Doppler effect** over a catalog you should use the function: 
```
doppler_boost(beta_var, beta_lat, beta_long, vecs, z)
```
>beta_var $\rightarrow$ absolute value of beta
>
>beta_lat $\rightarrow$ latitude of beta
>
>beta_long $\rightarrow$ longitude of beta
>
>vecs $\rightarrow$ the array with the x, y and z components of the direction vector for each object 
>
>z $\rightarrow$ is the observed redshifts' array for each object

---
