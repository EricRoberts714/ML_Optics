import numpy
import os

# input_array = numpy.loadtxt("https://raw.githubusercontent.com/EricRoberts714/ML_Optics/master/oasys_scripts/aps_axo_influence_functions2019.dat")
# basis = numpy.loadtxt("https://raw.githubusercontent.com/EricRoberts714/ML_Optics/master/oasys_scripts/aps_axo_orthonormal_functions2019.dat")

from calculate_simulated_deformations import write_surface_file

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
                       correction_file_M3=None,
                       calculate_correction=True):

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

    #
    # M3 corrector (not needed)
    #

    if calculate_correction:
        input_wavefront = output_wavefront


        output_wavefront, target_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_corrector1D(
            input_wavefront, grazing_angle=0.02181, focus_at=2.64, apodization=0, apodization_ratio=6.0,
            write_correction_profile=1)


    # from srxraylib.plot.gol import plot
    # plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity())


    if correction_file_M3 is not None:
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


if __name__ == "__main__":

    from srxraylib.plot.gol import plot, plot_table
    import matplotlib.pylab as plt

    import h5py

    file = h5py.File("set1_profiles.h5", 'r')
    subgroup_name = "profiles_stack"
    image_x = file[subgroup_name + "/X"][()]
    x0 = file[subgroup_name + "/Y"][()]
    PROFILES = file[subgroup_name + "/Z"][()].T

    print(PROFILES.shape, image_x.shape, x0.shape)

    correction_number_of_points = 200
    CORRECTIONS = numpy.zeros((image_x.size, correction_number_of_points))
    do_plot = False

    for entry_number in range(image_x.size):
        print(">>>>>>>>>>>>>>>>>>>>>>>>> running entry index: %d " % entry_number)

        error_file_M1 = "deformation.dat"
        f = open(error_file_M1, 'w')
        for i in range(x0.size):
            f.write("%g %g\n" % (x0[i], PROFILES[entry_number,i]))
        f.close()

        # make correction
        output_wavefront_zero = run_whole_beamline(error_file_M1=error_file_M1,
                           correction_file_M3=None, calculate_correction=False)

        output_wavefront_ideal = run_whole_beamline(error_file_M1=error_file_M1,
                           correction_file_M3=None, calculate_correction=True)

        output_wavefront_file = run_whole_beamline(error_file_M1=error_file_M1,
                           correction_file_M3="correction_profile1D.dat", calculate_correction=False)

        # interpolate to reduce the number of points to 200
        b = numpy.loadtxt("correction_profile1D.dat")
        x = numpy.linspace(b[0, 0], b[-1, 0], correction_number_of_points)
        y = numpy.interp(x, b[:, 0], b[:, 1], )

        CORRECTIONS[entry_number, :] = y

        if do_plot:
            plot(output_wavefront_zero.get_abscissas(), output_wavefront_zero.get_intensity(),
                 output_wavefront_file.get_abscissas(), output_wavefront_file.get_intensity(),
                 legend=["zero correction","perfectly corrected (from file)"],
                 title=error_file_M1,
                 show=0)
            plt.savefig(os.path.splitext(error_file_M1)[0]+".png")
            plt.show()

        if entry_number == 0:
            IMAGE_UNCORRECTED = numpy.zeros((image_x.size, output_wavefront_zero.get_abscissas().size ))
            IMAGE_CORRECTED = numpy.zeros((image_x.size, output_wavefront_file.get_abscissas().size))

        IMAGE_UNCORRECTED[entry_number, :] = output_wavefront_zero.get_intensity()
        IMAGE_CORRECTED[entry_number, :] = output_wavefront_file.get_intensity()


    write_surface_file(CORRECTIONS.T, image_x, x, "set1_corrections.h5", overwrite=True)
    write_surface_file(IMAGE_UNCORRECTED.T, image_x, output_wavefront_zero.get_abscissas(),
                       "set1_image_uncorrected.h5", overwrite=True)
    write_surface_file(IMAGE_CORRECTED.T, image_x,  output_wavefront_file.get_abscissas(),
                       "set1_image_corrected.h5", overwrite=True)
