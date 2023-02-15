import numpy as np
import h5py as h5
from functions import slope_fit

input_file_path = "/home/enda/PhD/Shell_Model/Data/ThesisData/SIM_DATA_[HYDRO-INTFACRK4-FULL]_N[25]_T[0.0,0.0001,1000000.000]_NU[5e-07]_ALPHA[1.500]_K[0.050,2.000]_EPS[0.50]_FORC[DELTA,1,0.100]_u0[N_SCALING]_TAG[Update]/"

with h5.File(input_file_path + "Stats_HDF_Data.h5", 'r') as in_file:
    print(list(in_file.keys()))
    str_func_vel_flux = in_file["StructureFunctionVelFluxAbs"][:, :, :]
    str_func_enrg_flux = str_func_vel_flux[:, :, 0]

with h5.File(input_file_path + "System_Measure_HDF_Data.h5", 'r') as in_file:
    print(list(in_file.keys()))
    k = in_file["k"][:]
    a_n_t_avg = in_file["TimeAveragedVelocityAmplitudes"][:]
    enrg_spec_t_avg = in_file["TimeAveragedEnergySpectrum"][:]


shell_start = np.arange(0, 10)
length = np.arange(6, 15)

# slope_data_list = []
# resid = []
min_slope_err = 1e10
min_resid = 1e10
for i in shell_start:
    for j in length:
        print("Shells: [{}, {}] - len: {} - ".format(i, i + j, j), end=" ")
        slope, c, resid = slope_fit(np.log(k), np.log(str_func_enrg_flux[:, 2]), i, i + j)
        print("slope: {} -- res: {}".format(slope, resid))
        slope_err = np.absolute(np.absolute(slope) - 1.000)
        if slope_err < min_slope_err:
            min_slope_err = slope_err
            slope_data_list = [i, i + j, np.absolute(slope), resid, np.absolute(np.absolute(slope) - 1.000)]
        if resid < min_resid:
            min_resid = resid
            resid_data_list = [i, i + j, np.absolute(slope), resid, np.absolute(np.absolute(slope) - 1.000)]

print(slope_data_list)
print(resid_data_list)


slope, c, resid = slope_fit(np.log(k), np.log(a_n_t_avg), slope_data_list[0], slope_data_list[1])
print(slope, resid)

slope, c, resid = slope_fit(np.log(k), np.log(enrg_spec_t_avg), slope_data_list[0], slope_data_list[1])
print(slope, resid)

slope, c, resid = slope_fit(np.log(k), np.log(a_n_t_avg), resid_data_list[0], resid_data_list[1])
print(slope, resid)