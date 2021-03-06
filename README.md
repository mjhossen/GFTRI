# GFTRI
GFTRI is a tsunami source inversion method based on time reverse imaging. It stands for "Green's Function based time Reverse Imaging". This method is defferent than the traditional least squares method. Instead of solving the evolved linear system by least squares method, it uses Green's functions (GFs) to estimate source amplitude by convolving them with observed waveforms in reversed time order.  For details, see the following articles. The current version of the code has been used for estimating sea surface displacement due to the 2018 Kodiak earthquake tsunami. It requires fault parameters, GFs and observation data stored in .h5 file. The data files are put in the Input folder except GFs file. GFs file cannot be uploaded to the repository due to it's large size. If you want to run the code with our data set, please contact me through the email: md.hossen@colorado.edu or mjhossen55@gmail.com. 


Hossen, M. J., Sheehan, A. F., & Satake, K. (2020). A multi-fault model estimation from tsunami data: An application to the 2018 M7. 9 Kodiak earthquake. Pure and Applied Geophysics, 1-12.

Hossen, M. J., Cummins, P. R., Dettmer, J., & Baba, T. (2015). Time reverse imaging for far‐field tsunami forecasting: 2011 Tohoku earthquake case study. Geophysical Research Letters, 42(22), 9906-9915.
