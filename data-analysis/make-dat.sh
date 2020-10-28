#!/usr/bin/env bash

#
# Preprocess data used in flowtime fit and global fit
#



echo "Performing Fourier transform on 2-point functions..."
find ./data/raw_data -name *wilson_twopt*.h5 -execdir latan-sample-ft -o ../../../../processed_data/fft_correlator/fft_{} {} \; 1>/dev/null

echo "Preprocessing c3 data against flowtime..."
for L in 64 128 256; do
emt-write-flowtime-c3 0.1 ${L} -0.0305
 /Users/s1854162/Desktop/wilsonflow/code_package/scalar_emt_analysis/src/build/emt-write-flowtime-c3 0.1 ${L} -0.031
 /Users/s1854162/Desktop/wilsonflow/code_package/scalar_emt_analysis/src/build/emt-write-flowtime-c3 0.2 ${L} -0.061
 /Users/s1854162/Desktop/wilsonflow/code_package/scalar_emt_analysis/src/build/emt-write-flowtime-c3 0.2 ${L} -0.062
 /Users/s1854162/Desktop/wilsonflow/code_package/scalar_emt_analysis/src/build/emt-write-flowtime-c3 0.3 ${L} -0.091
 /Users/s1854162/Desktop/wilsonflow/code_package/scalar_emt_analysis/src/build/emt-write-flowtime-c3 0.3 ${L} -0.092
done
  
echo "Generating renormalised mass data..."
latan-sample-fake  0.08408 0.00038 2000 data/processed_data/global_fit_data/mass/g0.1_m2-0.0305_mass_sample.h5
latan-sample-fake  0.03408 0.00038 2000 data/processed_data/global_fit_data/mass/g0.1_m2-0.031_mass_sample.h5 
latan-sample-fake  0.032435 0.000245 2000 data/processed_data/global_fit_data/mass/g0.2_m2-0.061_mass_sample.h5
latan-sample-fake  0.007435 0.000245 2000 data/processed_data/global_fit_data/mass/g0.2_m2-0.062_mass_sample.h5
latan-sample-fake  0.021505556 0.000173333 2000 data/processed_data/global_fit_data/mass/g0.3_m2-0.091_mass_sample.h5
latan-sample-fake  0.010394444 0.000173333 2000 data/processed_data/global_fit_data/mass/g0.3_m2-0.092_mass_sample.h5
