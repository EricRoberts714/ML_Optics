import numpy


def calculate_wavefront1D(wavelength=1e-10,
                          undulator_length=1.0, undulator_distance=10.0,
                          x_min=-0.1, x_max=0.1, number_of_points=101,
                          wavefront_position=0, add_random_phase=0):
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

    sigma_r = 2.740 / 4 / numpy.pi * numpy.sqrt(wavelength * undulator_length)
    sigma_r_prime = 0.69 * numpy.sqrt(wavelength / undulator_length)

    wavefront1D = GenericWavefront1D.initialize_wavefront_from_range(x_min=x_min, x_max=x_max,
                                                                     number_of_points=number_of_points)
    wavefront1D.set_wavelength(wavelength)

    if wavefront_position == 0:  # Gaussian source
        wavefront1D.set_gaussian(sigma_x=sigma_r, amplitude=1.0, shift=0.0)
    elif wavefront_position == 1:  # Spherical source, Gaussian intensity
        wavefront1D.set_spherical_wave(radius=undulator_distance, center=0.0, complex_amplitude=complex(1, 0))
        # weight with Gaussian
        X = wavefront1D.get_abscissas()
        A = wavefront1D.get_complex_amplitude()
        sigma = undulator_distance * sigma_r_prime
        sigma_amplitude = sigma * numpy.sqrt(2)
        Gx = numpy.exp(-X * X / 2 / sigma_amplitude ** 2)
        wavefront1D.set_complex_amplitude(A * Gx)

    if add_random_phase:
        wavefront1D.add_phase_shifts(2 * numpy.pi * numpy.random.random(wavefront1D.size()))

    return wavefront1D


def calculate_output_wavefront_after_reflector1D(input_wavefront, radius=10000.0, grazing_angle=1.5e-3, error_flag=0,
                                                 error_file="", write_profile=0):
    import numpy
    from scipy import interpolate
    output_wavefront = input_wavefront.duplicate()
    abscissas = output_wavefront.get_abscissas()
    abscissas_on_mirror = abscissas / numpy.sin(grazing_angle)
    if radius >= 0:
        height = radius - numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)
    else:
        height = radius + numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)

    if error_flag:
        a = numpy.loadtxt(error_file)
        finterpolate = interpolate.interp1d(a[:, 0], a[:, 1], fill_value="extrapolate")
        height_interpolated = finterpolate(abscissas_on_mirror)
        height += height_interpolated

    phi = -2 * output_wavefront.get_wavenumber() * height * numpy.sin(grazing_angle)
    output_wavefront.add_phase_shifts(phi)

    # output files
    if write_profile:
        f = open("reflector_profile1D.dat", "w")
        for i in range(height.size):
            f.write("%g %g\n" % (abscissas_on_mirror[i], height[i]))
        f.close()
        print("File reflector_profile1D.dat written to disk.")

    return output_wavefront, abscissas_on_mirror, height

def calculate_output_wavefront_after_corrector1D(input_wavefront, grazing_angle=1.5e-3, focus_at=10.0,
                                                 apodization=0, apodization_ratio=0.1, write_correction_profile=0):
    import numpy
    from scipy import interpolate
    output_wavefront = input_wavefront.duplicate()
    target_wavefront = input_wavefront.duplicate()
    target_wavefront.set_spherical_wave(radius=-focus_at, center=0.0, complex_amplitude=1.0)

    phase_input = input_wavefront.get_phase(unwrap=True)
    phase_target = target_wavefront.get_phase(unwrap=True)
    phase_correction = phase_target - phase_input
    abscissas = target_wavefront.get_abscissas()
    abscissas_on_mirror = abscissas / numpy.sin(grazing_angle)

    # output_wavefront.add_phase_shift(phase_correction)
    height = - phase_correction / (2 * output_wavefront.get_wavenumber() * numpy.sin(grazing_angle))

    if apodization == 0:
        height -= height[height.size // 2]
    elif apodization == 1:
        apodization = input_wavefront.get_intensity()
        apodization = (apodization / apodization.max())
        height *= apodization
        height -= height[0]
    elif apodization == 2:
        sigma = numpy.abs(abscissas[-1] - abscissas[0]) * apodization_ratio
        apodization = numpy.exp(- abscissas ** 2 / 2 / sigma ** 2)
        apodization /= apodization.max()
        height *= apodization
        height -= height[0]

    # calculate phase shift from new profile
    phi = -2 * output_wavefront.get_wavenumber() * height * numpy.sin(grazing_angle)
    output_wavefront.add_phase_shift(phi)

    # output files
    if write_correction_profile:
        f = open("correction_profile1D.dat", "w")
        for i in range(height.size):
            f.write("%g %g\n" % (abscissas_on_mirror[i], height[i]))
        f.close()
        print("File correction_profile1D.dat written to disk.")

    return output_wavefront, target_wavefront, abscissas_on_mirror, height



def run_whole_beamline(error_file_M1="/home/manuel/Oasys/dabam_profile_140327232022424.dat",
                       correction_file_M3="/home/manuel/Oasys/dabam_profile_140522663056440.dat"):
    #
    # source
    #
    output_wavefront = calculate_wavefront1D(wavelength=4.9593679373280105e-09,
                                             wavefront_position=1,
                                             undulator_length=3.98,
                                             undulator_distance=13.73,
                                             x_min=-0.00147,
                                             x_max=0.00147,
                                             number_of_points=3000,
                                             add_random_phase=False)





    #
    # M1
    #
    input_wavefront = output_wavefront
    output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_reflector1D(input_wavefront,radius=496600.0,grazing_angle=0.0218,error_flag=1,error_file=error_file_M1,write_profile=1)



    #
    # screen at M3 position
    #

    #
    # Import section
    #
    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D


    input_wavefront = output_wavefront

    #
    # info on current oe
    #
    #
    #    -------WOScreen1D---------
    #        -------BoundaryShape---------
    #

    #
    # define current oe
    #
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    #
    # propagating (***  ONLY THE ZOOM PROPAGATOR IS IMPLEMENTED ***)
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=0.000000, q=13.599000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),
                                                   propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 2.000000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')




    #
    # M3 circlular mirror
    #



    input_wavefront = output_wavefront
    output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_reflector1D(input_wavefront,
                                                                                                 radius=220.71532,
                                                                                                 grazing_angle=0.02181,
                                                                                                 error_flag=0,
                                                                                                 error_file="",
                                                                                                 write_profile=0)

    # #
    # # M3 corrector (not needed)
    # #
    #
    # input_wavefront = output_wavefront
    #
    #
    # output_wavefront, target_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_corrector1D(
    #     input_wavefront, grazing_angle=0.02181, focus_at=2.64, apodization=0, apodization_ratio=6.0,
    #     write_correction_profile=1)
    #
    #
    # from srxraylib.plot.gol import plot
    # plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity())


    #
    # M3 corrector from external file
    #





    input_wavefront = output_wavefront
    output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_reflector1D(input_wavefront,
                                                                                                 radius=100000000.0,
                                                                                                 grazing_angle=0.02181,
                                                                                                 error_flag=1,
                                                                                                 error_file=correction_file_M3,
                                                                                                 write_profile=0)

    #
    # propagation to sample
    #



    input_wavefront = output_wavefront


    #
    # define current oe
    #
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    #
    # propagating (***  ONLY THE ZOOM PROPAGATOR IS IMPLEMENTED ***)
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=0.000000, q=2.640000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),
                                                   propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.020000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    return output_wavefront

def get_surface_from_basis_and_coefficients(basis,coefficients= \
        [-1.34111718e-08 ,2.89091778e-11 , 4.55787593e-07,  3.71667246e-08,
          1.03623232e-07,  3.30400686e-08, -3.41480067e-09,  5.67368681e-10,
         -1.08934818e-08, -1.67860465e-08, -9.94266810e-09, -1.41912905e-08,
         -1.32628413e-08, -1.36259836e-08, -9.41153403e-09, -9.47141468e-09,
         -5.08842152e-09, -6.38382702e-10, -5.57004604e-09, -2.95162052e-09,],
                                            ):
    out = 0
    if isinstance(coefficients,list):
        coefficients = numpy.array(coefficients)

    for i in range(coefficients.size):
        out += coefficients[i] * basis[:,i]
    return out




if __name__ == "__main__":




    from srxraylib.plot.gol import plot


    #
    # run using the perfect correction profile
    #

    output_wavefront = run_whole_beamline(error_file_M1="deformation.dat",
                       correction_file_M3="correction.dat")

    plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(),title="prefectly corrected")

    #
    # run using zero correction
    #


    x = numpy.linspace(-0.135,0.135,100)
    height = x * 0

    f = open("tmp.dat",'w')
    for i in range(x.size):
        f.write("%g %g \n"%(x[i],height[i]))
    f.close()


    output_wavefront = run_whole_beamline(error_file_M1="deformation.dat",
                       correction_file_M3="tmp.dat")

    plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title="correction file tmp.dat: zero correction")



    #
    # run a correction from coefficients
    #




    # loads file with data to fit
    input_array = numpy.loadtxt("aps_axo_influence_functions2019.dat")
    abscissas = input_array[:, 0].copy()
    basis = numpy.loadtxt("aps_axo_orthonormal_functions2019.dat")
    #
    print(">>>>basis",basis.shape)
    coefficients = [-1.34111718e-08, 2.89091778e-11, 4.55787593e-07, 3.71667246e-08, \
     1.03623232e-07, 3.30400686e-08, -3.41480067e-09, 5.67368681e-10, \
     -1.08934818e-08, -1.67860465e-08, -9.94266810e-09, -1.41912905e-08, \
     -1.32628413e-08, -1.36259836e-08, -9.41153403e-09, -9.47141468e-09, \
     -5.08842152e-09, -6.38382702e-10, -5.57004604e-09, -2.95162052e-09, ]



    # coefficients = numpy.random.random(20) * 1e-9
    height0 = get_surface_from_basis_and_coefficients(basis,coefficients)



    # abscissas for the mirror
    x = numpy.linspace(-0.135,0.135,100)
    height = numpy.interp(x,abscissas/1000,height0)
    print(x)
    #
    plot(abscissas/1000,height0,x,height*1.1, legend=["evaluated from coeffcients","interpolated to the mirror coordinates"])
    # plot(abscissas/1000, height0)
    # plot(x, height)

    f = open("tmp.dat",'w')
    for i in range(x.size):
        f.write("%g %g \n"%(x[i],height[i]))
    f.close()


    output_wavefront = run_whole_beamline(error_file_M1="deformation.dat",
                       correction_file_M3="tmp.dat")

    plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title="correction file tmp.dat: from coefficients")





