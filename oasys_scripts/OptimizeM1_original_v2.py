

import numpy
from flexon_ken_memo_optimization import run_whole_beamline, get_surface_from_basis_and_coefficients
from srxraylib.plot.gol import plot

from scipy import optimize
import nevergrad as ng
import time

## Define function to minimize/maximize ##
##########################################

def get_corrector_profile(coefficients):
    # This function takes as input the 1 by 20 list of M1 correction coefficients
    # and returns the maximum intensity of the beamline corrected by the
    # coefficients. The negative of the intensity is returned so that the
    # minimization algorithms chosen will optimize with respect to the max
    # intensity.

    ## Load data and basis
    #input_array = numpy.loadtxt("aps_axo_influence_functions2019.dat")
    abscissas = input_array[:, 0].copy()
    #basis = numpy.loadtxt("aps_axo_orthonormal_functions2019.dat")

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


    ## Get max intensity ##
    #######################
    intensity = output_wavefront.get_intensity()
    abscissas = output_wavefront.get_abscissas()

    maxIndex = numpy.argmax(intensity)
    maxIntensity = intensity[maxIndex]
    distance = numpy.abs( abscissas[maxIndex] ) * 1e6
    #maxIntensityTrue = 212.97657172329227  # The optimally corrected max
                                            # intensity from earlier testing

    #print(maxIntensityTrue - maxIntensity)
    # return(-maxIntensity)
    return(-maxIntensity / (1 + distance))


if __name__ == "__main__":
    ## Original Scripts (Perfectly Corrected) ##
    ############################################

    output_wavefront0 = run_whole_beamline(error_file_M1="deformation.dat",
                                          correction_file_M3="correction.dat")


    input_array = numpy.loadtxt("aps_axo_influence_functions2019.dat")
    abscissas = input_array[:, 0].copy()
    basis = numpy.loadtxt("aps_axo_orthonormal_functions2019.dat")

    print(abscissas)

    aaa = numpy.loadtxt('correction.dat')
    x00 = aaa[:,0]
    height00 = aaa[:,1]
    # plot(x00,height00*1.0, legend=["Perfect corrector"])

    ## Final beam
    #plot(output_wavefront0.get_abscissas(), output_wavefront0.get_intensity(), title="correction file tmp.dat: from coefficients")


    coefficients = [-1.34111718e-08, 2.89091778e-11, 4.55787593e-07, 3.71667246e-08,  1.03623232e-07, 3.30400686e-08, -3.41480067e-09, 5.67368681e-10,  -1.08934818e-08, -1.67860465e-08, -9.94266810e-09, -1.41912905e-08,  -1.32628413e-08, -1.36259836e-08, -9.41153403e-09, -9.47141468e-09,  -5.08842152e-09, -6.38382702e-10, -5.57004604e-09, -2.95162052e-09, ]

    height0 = get_surface_from_basis_and_coefficients(basis,coefficients)
    ## Interpolate here
    x = numpy.linspace(-0.135,0.135,100)
    height = numpy.interp(x,abscissas/1000,height0)
    # plot(x00,height00*1.0,
    #      abscissas/1000,height0,
    #      x,height*1.0,
    #      legend=["Perfect corrector","evaluated from coeffcients","interpolated to the mirror coordinates"])

    print((x.size))

    f = open("tmp.dat",'w')
    for i in range(x.size):
        f.write("%g %g \n"%(x[i],height[i]))
    f.close()


    output_wavefront = run_whole_beamline(error_file_M1="deformation.dat",
                       correction_file_M3="tmp.dat")

    # plot(output_wavefront0.get_abscissas(), output_wavefront0.get_intensity(),
    #     output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
    #      title="correction file tmp.dat: from coefficients",legend=["theoretical","from coefs"])

    maxVal = numpy.amax(output_wavefront.get_intensity())
    maxValTrue = 212.97657172329227             # From earlier testing
    maxValNoCorrection = 111.43902895363388     # From earlier testing
    print('The max value is ', maxVal)
    print('The true max val ', maxValTrue)
    print('Max, no correction ', maxValNoCorrection)




    start_time = time.time()

    # method = 0  # good
    # method = 1  # fail
    # method = 2  # good

    method = 10  # nelder-mead good
    # method = 11  # L-BFGS-B fail
    # method = 12 # powell bad
    # method = 13 # basinhopping fail


    if method == 0:
        optimizer = ng.optimizers.OnePlusOne(parametrization = 20, budget = 1000)
        recommendation = optimizer.minimize(get_corrector_profile)
        print(recommendation.value)
        coefficientsOptimized = recommendation.value
    elif method == 1:
        optimizer = ng.optimizers.CM(parametrization = 20, budget = 5000)
        recommendation = optimizer.minimize(get_corrector_profile)
        print(recommendation.value)
        coefficientsOptimized = recommendation.value
    elif method == 2:
        optimizer = ng.optimizers.OnePlusOne(parametrization = 20, budget = 500)
        recommendation = optimizer.minimize(get_corrector_profile)
        print(recommendation.value)
        coefficientsOptimized = recommendation.value
    elif method == 10:
        coefficients = numpy.zeros(20)
        result = optimize.minimize(get_corrector_profile, coefficients, method='nelder-mead',
                                   options={'maxiter': 5000, 'fatol': 1e-4, 'xatol': 1e-4, 'disp': True})
        coefficientsOptimized = result.x
    elif method == 11:
        coefficients = numpy.zeros(20)
        result = optimize.minimize(get_corrector_profile, coefficients, method='L-BFGS-B',
                                   options={'maxiter': 5000, 'ftol': 1e-4, 'gtol': 1e-4, 'disp': True})
        coefficientsOptimized = result.x
    elif method == 12:
        coefficients = numpy.zeros(20)
        result = optimize.minimize(get_corrector_profile, coefficients, method='Powell',
                                   options={'maxiter': 5000, 'ftol': 1e-4, 'xtol': 1e-4, 'disp': True})
        coefficientsOptimized = result.x
    elif method == 13:
        coefficients = numpy.zeros(20)
        result = optimize.basinhopping(get_corrector_profile, coefficients, niter=10)
        coefficientsOptimized = result.x

    NGtime = time.time() - start_time
    print(NGtime)


    height0 = get_surface_from_basis_and_coefficients(basis, coefficientsOptimized)

    import matplotlib.pylab as plt

    plt.bar(numpy.arange(20)-0.1, numpy.array(coefficients), width=0.2)
    plt.bar(numpy.arange(20) + 0.1, numpy.array(coefficientsOptimized), width=0.2)
    plt.title("Coeff OK - Optimized")
    plt.show()

    ## Interpolate here
    x = numpy.linspace(-0.135,0.135,100)
    height = numpy.interp(x,abscissas/1000,height0)

    plot(x00,height00*1.0,
         abscissas/1000,height0,
         x,height*1.0,
         legend=["perfect corrector","evaluated from coeffcients","interpolated to the mirror coordinates"])

    f = open("tmp.dat",'w')
    for i in range(x.size):
        f.write("%g %g \n"%(x[i],height[i]))
    f.close()

    correction_file_M3 = "tmp.dat"

    output_wavefront = run_whole_beamline(error_file_M1="deformation.dat",
                       correction_file_M3="tmp.dat")

    maxVal = numpy.amax(output_wavefront.get_intensity())
    print('Max intensity after minimizing is:', maxVal)

    plot(output_wavefront0.get_abscissas(), output_wavefront0.get_intensity(),
         output_wavefront.get_abscissas(), output_wavefront.get_intensity(),
         legend=["Theoretical","From minimization"],
         title="correction file tmp.dat: from coefficients")
