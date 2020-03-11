## Imports ##
#############

import numpy
from flexon_ken_memo_optimization import run_whole_beamline, get_surface_from_basis_and_coefficients
from srxraylib.plot.gol import plot

from scipy import optimize
import time

## Define function to minimize/maximize ##
##########################################
def get_corrector_profile(coefficients,do_plot=0,return_wavefront=False):
    # This function takes as input the 1 by 20 list of M1 correction coefficients
    # and returns the maximum intensity of the beamline corrected by the
    # coefficients. The negative of the intensity is returned so that the
    # minimization algorithms chosen will optimize with respect to the max
    # intensity.

    ## Load data and basis
    input_array = numpy.loadtxt("aps_axo_influence_functions2019.dat")
    abscissas = input_array[:, 0].copy()
    basis = numpy.loadtxt("aps_axo_orthonormal_functions2019.dat")

    ## Get height values and interpolate grid
    height0 = get_surface_from_basis_and_coefficients(basis, coefficients)
    x = numpy.linspace(-0.135, 0.135, 100)
    height = numpy.interp(x, abscissas / 1000, height0)
    #print(">>>>>", height.size)
    # plot(abscissas/1000,height0,x,height*1.0, legend=["evaluated from coeffcients","interpolated to the mirror coordinates"])

    ## Write data in tmp file and assign to correction file
    f = open("tmp.dat", 'w')
    for i in range(x.size):
        f.write("%g %g \n" % (x[i], height[i]))
    f.close()
    correction_file_M3 = "tmp.dat"

    ## Run correcton
    output_wavefront = run_whole_beamline(error_file_M1="deformation.dat",
                                          correction_file_M3="tmp.dat")

    if do_plot:
        from srxraylib.plot.gol import plot
        plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title=str(coefficients))
    ## Get max intensity ##
    #######################
    maxIntensity = numpy.amax(output_wavefront.get_intensity())
    #maxIntensityTrue = 212.97657172329227  # The optimally corrected max
                                            # intensity from earlier testing

    #print(maxIntensityTrue - maxIntensity)
    if return_wavefront:
        return output_wavefront
    else:
        return(-maxIntensity)

## Define coefficients ##
#########################



numpy.random.seed(12345)
coefficientsRandom = numpy.random.random(20) * 1e-7     #Random sampling

wf_ramdom = get_corrector_profile(coefficientsRandom,do_plot=0,return_wavefront=1)

coefficients = [-1.34111718e-08, 2.89091778e-11, 4.55787593e-07, 3.71667246e-08, \
    1.03623232e-07, 3.30400686e-08, -3.41480067e-09, 5.67368681e-10, \
    -1.08934818e-08, -1.67860465e-08, -9.94266810e-09, -1.41912905e-08, \
    -1.32628413e-08, -1.36259836e-08, -9.41153403e-09, -9.47141468e-09, \
    -5.08842152e-09, -6.38382702e-10, -5.57004604e-09, -2.95162052e-09, ]

wf_ideal = get_corrector_profile(coefficients,do_plot=0,return_wavefront=1)





# tmp = get_surface_from_coefficients(coefficients,do_plot=1,title="ideal",show=1)
## Optimize ##
##############

# We use a simple gradient descent method title 'Nelder-Mead' which requires no Jacobian
# or Hessian

# Define minimization scheme parameters
maxiter = 10000
fatol = 1e-3    # Absolute error in input variables x between iterations
xatol = 1e-3    # Absolute error in f(x) between iterations

# Random coefficient result
resultRandom = optimize.minimize(get_corrector_profile, coefficientsRandom, method='nelder-mead',options = {'maxiter':maxiter,'fatol':fatol, 'xatol':xatol, 'disp':True})
print(resultRandom)     #or print(resultRandom.x)

wf_optimized = get_corrector_profile(resultRandom["x"],do_plot=1,return_wavefront=1)

plot(wf_ramdom.get_abscissas(),wf_ramdom.get_intensity(),
     wf_ideal.get_abscissas(),wf_ideal.get_intensity(),
     wf_optimized.get_abscissas(),wf_optimized.get_intensity(),
     legend=["random start","ideal","optimized"])


# # Sample coefficient result
# result = optimize.minimize(get_corrector_profile, coefficients, method='nelder-mead',options = {'maxiter':maxiter,'fatol':fatol, 'xatol':xatol, 'disp':True})
# print(result)     #or print(result.x)
# tmp = get_corrector_profile(result,do_plot=1)




