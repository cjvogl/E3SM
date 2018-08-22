from plot_convergence_test import plot_convergence_test

testName = 'dcmip2012_test41'
suffix = '_dcmip4_X1'
dtRef = 10.0
indRef = 15
varRef = 'u'
fileRef = './output_tsteptype5_tstep%2.0f_hydrostatic_dcmip4_X1/%s.nc' \
                  % (dtRef, testName)
plot_convergence_test(testName, fileRef, dtRef, indRef, varRef, suffix=suffix)
