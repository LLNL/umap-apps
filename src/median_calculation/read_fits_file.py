import astropy.io.fits
import sys
import math

# Reference
# http://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html#sphx-glr-generated-examples-io-plot-fits-image-py

# Usage
# python3 src/median_calculation/read_fits_file.py [x-slope] [x-intercept] [y-slope] [y-intercept] [#of input filts files] [prefix of the fits files] [path to a timestamp file]

def read_pixel(fits, x, y):
    #print (fits)
    image_data = astropy.io.fits.getdata(fits, ext=0)
    if (y >= 0 and y < image_data.shape[0]) and (x >= 0 and x < image_data.shape[1]):
        print(image_data[y][x])
    else:
        print(-1)


def read_timestamp(fname):
    time_stamp = []
    with open(fname) as f:
        for line in f:
            time_stamp.append(float(line))
    return time_stamp


def main(argv):
    x_slope = float(argv[1]);
    x_intercept = float(argv[2]);
    y_slope = float(argv[3]);
    y_intercept = float(argv[4]);

    fits_list = []
    for i in range(int(argv[5])):
        fits_list.append(argv[6] + str(i + 1) + ".fits")

    time_stamp = read_timestamp(argv[7])

    for k in range(len(fits_list)):
        offset = time_stamp[k] - time_stamp[0]
        x = int(round(x_slope * offset + x_intercept, 0))
        y = int(round(y_slope * offset + y_intercept, 0))
        print(x, ", ", y, ", ", k, end=' : ')
        read_pixel(fits_list[k], x, y)


if __name__ == '__main__':
    main(sys.argv)
