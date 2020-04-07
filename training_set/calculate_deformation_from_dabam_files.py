import numpy
from srxraylib.metrology.dabam import dabam
from srxraylib.plot.gol import plot


dm = dabam()

def load_dabam_profile(entry_number,
                       mirror_length = 2 * 0.07,
                       mirror_points = 100,
                       mirror_rms = 25e-9,
                       do_plot = True
                       ):

    # print(dm.inputs)

    dm.inputs['entryNumber'] = entry_number
    dm.load()  # access data

    # dm.inputs['plot'] = "heights"
    # dm.plot()

    x00 = dm.y.copy()
    y00 = dm.zHeights.copy()

    #center
    x0 = x00 - x00[x00.size//2]
    y0 = y00.copy()
    #rescale abscissas
    x0 = x0 / numpy.abs(x0).max() * 0.5 * mirror_length
    # rescale heights
    print("RMS mirror: %5.3f nm "%(1e9*y0.std()))
    y0 = y0 * mirror_rms / y0.std()
    # interpolate
    x = numpy.linspace(-0.5 * mirror_length, 0.5 * mirror_length, mirror_points)
    y = numpy.interp(x, x0, y0)
    if do_plot:
        plot(x00, y00,
             x0, y0,
             x, y,
             legend=["Original","Transformed","Interpolated"])
    return x,y

if __name__ == "__main__":

    for entry_number in range(1,83):
        print(">>>> entry: ",entry_number)
        x, y = load_dabam_profile(entry_number, do_plot=False)


        filename = "deformation%02d.dat" % entry_number
        f = open(filename, "w")
        for i in range(x.size):
            f.write("%g %g\n" % (x[i], y[i]))
        f.close()
        print("File written to disk: %s" % filename)

