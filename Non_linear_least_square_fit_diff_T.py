import numpy as np
from scipy.optimize import curve_fit
from scipy.special import lambertw
import os

# Define the parameters

d = 3  # dimension
np_f = 500000  # number of points for the fitted function
# np_d = 5  # number of data points
np_throw = [3, 7, 3, 7, 3, 3]
np_d = [6, 10, 6, 9, 6, 5]
#np_d_max = max(np_d)

output_dir='FITS'
try:
	os.mkdir(output_dir)
except:
	print("WARNING: output directory already exist")

def function_eps_star(eps_dot, alpha, tau0):
    temp1 = 2.0 * (d - 1.0)
    temp2 = temp1 * alpha * (eps_dot * tau0)**(-temp1)
    temp3 = temp1 * alpha
    obj = np.real(((lambertw(temp2)) / (temp3))**(-1.0 / temp1))

    temp4 = lambertw(temp2)

    return obj


file_pref = 'inputs/peak_pos_T'

#file_out_pref = 'FITS_Seperate_time_scale/'
file_out_pref = 'FITS/'
file_suf = ['04', '06', '065', '07', '075', '08']
Temperature = [0.40, 0.60, 0.65, 0.70, 0.75, 0.80]
N_particle = [3328, 3328, 3328, 3328, 3328, 3328]

# with open("Input", 'r') as f1:
for inum, n_d in enumerate(np_d):
    eps_dot = np.zeros((n_d), dtype=np.float64)
    eps_star_data = np.zeros((n_d), dtype=np.float64)
    file_name = file_pref + file_suf[inum] + '.dat'
    with open(file_name, 'r') as f1:
        for i in range(np_throw[inum]):
            f1_temp = f1.readline()

        for i in range(n_d):
            f1_temp = f1.readline()
            a, b = f1_temp.split()
            eps_dot[i] = float(a)  # * 10.0**15
            eps_star_data[i] = float(b)
#            print(eps_dot[i], eps_star_data[i])

    bnds = ((0, [np.inf, np.inf]))

    popt, pcov = curve_fit(function_eps_star, eps_dot,
                           eps_star_data, bounds=bnds)

    alpha_f = popt[0]
    tau0_f = popt[1]
    #pcov = params[1, :]
    perr = np.sqrt(np.diag(pcov))
    alpha_err = perr[0]
    tau0_err = perr[1]

##    print("perr", perr)

#eps_dot_i = np.zeros((10, 1), dtype=np.float64)
    eps_dot_i = np.logspace(-20, -1, np_f)


    file_out_1 = file_out_pref + 'FITTED_T_' + file_suf[inum]
    fi=open(file_out_1, "w")
    for i in range(np_f):
        eps_star = function_eps_star(eps_dot_i[i], alpha_f, tau0_f)
        print(eps_dot_i[i], eps_star, file=fi)
    fi.close()


    file_out_2 = file_out_pref + 'MASTER_T_' + file_suf[inum]
    f2=open(file_out_2, "w")

    file_out_3 = file_out_pref + 'T_alpha_tau0_errors_N'
    f3=open(file_out_3, "w")
    for i2 in range(n_d):
        print((tau0_f * eps_dot[i2] * (4.0 * alpha_f)
                   ** (-0.25)), (alpha_f / eps_star_data[i2]**4), eps_dot[i2], eps_star_data[i2], file=f2)

        print(Temperature[inum], alpha_f, tau0_f, alpha_err,
              tau0_err, N_particle[inum], file=f3)

    f2.close()
    f3.close()
