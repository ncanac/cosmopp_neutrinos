import os
import numpy as np
from classy import Class

class bao_boss:

    # initialization routine

    def __init__(self, datafile):
        # define array for values of z and data points
        self.z = np.array([], 'float64')
        self.data = np.array([], 'float64')
        self.error = np.array([], 'float64')
        self.type = np.array([], 'int')

        # read redshifts and data points
        with open(datafile, 'r') as filein:
            for line in filein:
                if line.strip() and line.find('#') == -1:
                    # the first entry of the line is the identifier
                    this_line = line.split()
                    self.z = np.append(self.z, float(this_line[1]))
                    self.data = np.append(self.data, float(this_line[2]))
                    self.error = np.append(self.error, float(this_line[3]))
                    self.type = np.append(self.type, int(this_line[4]))

        # number of data points
        self.num_points = np.shape(self.z)[0]

        # end of initialization

    # compute likelihood

    def loglkl(self, params):
        cosmo = Class()
        cosmo.set(params)
        cosmo.compute()

        chi2 = 0.

        # for each point, compute angular distance da, radial distance dr,
        # volume distance dv, sound horizon at baryon drag rs_d,
        # theoretical prediction and chi2 contribution
        for i in range(self.num_points):

            da = cosmo.angular_distance(self.z[i])
            dr = self.z[i] / cosmo.Hubble(self.z[i])
            dv = pow(da * da * (1 + self.z[i]) * (1 + self.z[i]) * dr, 1. / 3.)
            rs = cosmo.rs_drag()

            if self.type[i] == 3:
                theo = dv / rs

            elif self.type[i] == 4:
                theo = dv

            elif self.type[i] == 5:
                theo = da / rs

            elif self.type[i] == 6:
                theo = 1. / cosmo.Hubble(self.z[i]) / rs

            elif self.type[i] == 7:
                theo = rs / dv

            chi2 += ((theo - self.data[i]) / self.error[i]) ** 2

        # return ln(L)
        # lkl = - 0.5 * chi2
        # return -2ln(L)
        lkl = chi2

        return lkl
