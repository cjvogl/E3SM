from plot_convergence_test import plot_convergence_test

testName = 'dcmip2012_test31'
suffix = '_nu0.0'
dtRef = 0.000390625
nmax = 768000
indRef = 10
varRef = 'T'
fileRef = './output_tsteptype5_tstep%10.9f_nmax%d_nu0.0/%s.nc' \
              % (dtRef, nmax, testName)
plot_convergence_test(testName, fileRef, dtRef, indRef, varRef, \
                          minDt=0.07, maxDt=1.1, suffix=suffix)
