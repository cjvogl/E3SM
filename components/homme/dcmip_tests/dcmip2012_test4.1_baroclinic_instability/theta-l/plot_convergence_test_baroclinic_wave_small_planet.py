from plot_convergence_test import plot_convergence_test

testName = 'dcmip2012_test41'
suffix = '_dcmip4_X100'
dtRef = 0.1
indRef = 15
fileRef = './output_tsteptype5_tstep%2.1f_dcmip4_X100/%s.nc' \
              % (dtRef, testName)
plot_convergence_test(testName, fileRef, dtRef, indRef, varRef, suffix=suffix)
