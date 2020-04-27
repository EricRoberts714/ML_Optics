import numpy
from srxraylib.metrology.dabam import dabam
from srxraylib.plot.gol import plot, plot_image

import h5py
import os, time



def simulate_profile_1D_fractal(step=1.0, npoints=None, mirror_length=200.0,
                                power_law_exponent_beta=1.5,npoints_ratio_f_over_x=1.0, random_seed=8787,
                                renormalize_to_heights_sd=None,renormalize_to_slopes_sd=None,
                                frequency_max=None,frequency_min=None):
    """
    #
    # generates a 1-dimensional random rough surface z(x) with n_surface_points surface points.
    # The surface has a power lar PSD |f|**(-beta).
    # It is a fractal profile if 1<beta<3
    #
    :param step: step in mirror length (default=0.2)
    :param mirror_length: profile length
    :param npoints: number of points in mirror length (default=None, i.e., undefined so use step and mirror_length
                    to calculate it. If defined, use npoints and step is irrelevant)
    :param power_law_exponent_beta: beta value
    :param npoints_ratio_f_over_x: ratio of the number of points in frequency domain over real space (default=1.0)
    :param random_seed: a random seed to initialize numpy.seed(). Use zero to avoid initialization (default=8787)
    :param renormalize_to_heights_sd: set to a value to renormalize the profile to this height stdev value (default=None)
    :param renormalize_to_slopes_sd: set to a value to renormalize the profile to this slope stdev value (default=None)
    :param frequency_max:
    :param frequency_min:
    :return: (x,prof) where x = profile abscissas, prof = profile heights
    """

    if npoints is None:
        n_surface_points = int(1 + (mirror_length / step))
    else:
        n_surface_points = npoints

    if random_seed != 0:
        numpy.random.seed(seed=random_seed)


    x_coords = numpy.linspace(-0.5*mirror_length,0.5*mirror_length,n_surface_points)


    if frequency_min is None:
        f_from =  1/(1*mirror_length)
    else:
        f_from = frequency_min

    if frequency_max is None:
        f_to = 1/(2*step)
    else:
        f_to = frequency_max

    if npoints_ratio_f_over_x == 1.0:
        f_npoints = n_surface_points
    else:
        f_npoints = int(n_surface_points*npoints_ratio_f_over_x)

    freq = numpy.linspace(f_from,f_to,f_npoints)
    #todo: make exponent of power law a parameter
    ampl = freq**(-power_law_exponent_beta/2)
    phases = numpy.random.rand(freq.size)*2*numpy.pi
    ymirr = numpy.zeros(n_surface_points)
    for i in range(f_npoints):
        ymirr += (ampl[i] *  numpy.sin(2*numpy.pi*freq[i]*x_coords + phases[i]))

    if renormalize_to_heights_sd != None:
        ymirr = ymirr / ymirr.std() * renormalize_to_heights_sd

    if renormalize_to_slopes_sd != None:
        yslopes = numpy.gradient(ymirr, step)
        ymirr = ymirr / yslopes.std() * renormalize_to_slopes_sd

    return x_coords, ymirr

def write_surface_file(zz, xx, yy, file_name, overwrite=True, subgroup_name="profiles_stack"):

    if (os.path.isfile(file_name)) and (overwrite==True): os.remove(file_name)

    if not os.path.isfile(file_name):  # if file doesn't exist, create it.
        file = h5py.File(file_name, 'w')
        # points to the default data to be plotted
        file.attrs['default']          = subgroup_name + '/Z'
        # give the HDF5 root some more attributes
        file.attrs['file_name']        = file_name
        file.attrs['file_time']        = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        file.attrs['creator']          = 'write_surface_file'
        file.attrs['code']             = 'Oasys'
        file.attrs['HDF5_Version']     = h5py.version.hdf5_version
        file.attrs['h5py_version']     = h5py.version.version
        file.close()

    file = h5py.File(file_name, 'a')

    try:
        f1 = file.create_group(subgroup_name)
    except:
        f1 = file[subgroup_name]

    f1z = f1.create_dataset("Z", data=zz)
    f1x = f1.create_dataset("X", data=xx)
    f1y = f1.create_dataset("Y", data=yy)


    # NEXUS attributes for automatic plot
    f1.attrs['NX_class'] = 'NXdata'
    f1.attrs['signal'] = "Z"
    f1.attrs['axes'] = [b"Y", b"X"]

    f1z.attrs['interpretation'] = 'image'
    f1x.attrs['long_name'] = "X [m]"
    f1y.attrs['long_name'] = "Y [m]"


    file.close()

if __name__ == "__main__":

    number_of_profiles = 5000

    mirror_length = 2 * 0.07
    mirror_points = 100
    mirror_rms = 25e-9

    PROFILES = numpy.zeros((number_of_profiles, mirror_points))

    for i in range(number_of_profiles):

        power_law_exponent_beta = 1 + 2.0 * numpy.random.random() # in interval [1,3]

        x, y = simulate_profile_1D_fractal(step=mirror_length / (mirror_points - 1),
                                    npoints=None,
                                    mirror_length=mirror_length,
                                    power_law_exponent_beta=power_law_exponent_beta,
                                    npoints_ratio_f_over_x=1.0,
                                    random_seed=0,
                                    renormalize_to_heights_sd=mirror_rms,
                                    renormalize_to_slopes_sd=None,
                                    frequency_max=None,
                                    frequency_min=None)

        # print(x.shape, y.shape)
        PROFILES[i,:] = y
        # plot(x, y, title="beta=%f" % power_law_exponent_beta)

    image_x = numpy.arange(number_of_profiles)
    plot_image(PROFILES, image_x, x, aspect="auto")

    #
    # h5file
    #
    write_surface_file(PROFILES.T, image_x, x, "set1_profiles.h5", overwrite=True)
