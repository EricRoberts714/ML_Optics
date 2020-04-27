import numpy
from srxraylib.plot.gol import plot, plot_image
import h5py


def extract_profile(filename, index):
    file = h5py.File(filename, 'r')
    subgroup_name = "profiles_stack"
    image_x = file[subgroup_name + "/X"][()]
    x0 = file[subgroup_name + "/Y"][()]
    PROFILES = file[subgroup_name + "/Z"][()].T
    print(PROFILES.shape, image_x.shape, x0.shape)
    return x0, PROFILES[index, :].copy()

if __name__ == "__main__":

    profile_index = 111

    x, y = extract_profile("set1_profiles.h5", profile_index)
    xc, yc = extract_profile("set1_corrections.h5", profile_index)

    plot(x,  1e9 * y,
         xc, 1e9 * yc,
         xtitle="x [m]", ytitle="height [nm]", title="profile #%d" % profile_index,
         legend=["deformation in M1","correction in M3"])




    X, Y = extract_profile("set1_image_uncorrected.h5", profile_index)
    XC, YC = extract_profile("set1_image_corrected.h5", profile_index)

    plot(1e6 * X,  Y,
         1e6 * XC, YC,
         xtitle="x [um]", ytitle="intensity [a.u.]", title="profile #%d" % profile_index,
         legend=["uncorrected","corrected"])

    plot(x,  1e9 * y,
         xc / (xc[-1] - xc[0]) * (x[-1] - x[0]), -1e9 * yc,
         xtitle="x [m]", ytitle="height [nm]", title="profile #%d" % profile_index,
         legend=["deformation in M1","correction in M3 expanded and inverted"])