import yt
import numpy as np
import matplotlib.pyplot as plt
import sys

#sys.path.append('/home/jb5027/cfd-solver/')

# from src.mesh import Mesh
# from src.input import Input
# from src.mech_reader import MechReader
# from src.eos import temp_yk_to_e, temp_yk_to_h, rho_temp_yk_to_p
# from src.fluxes import get_fluxes
# from src.analysis import get_terms_for_analysis

# Load in two time steps
# ds_286092 = yt.load('/scratch/gpfs/jb5027/cavity-flame-data/C2/timesteps/plt286092')
# ds_286093 = yt.load('/scratch/gpfs/jb5027/cavity-flame-data/C2/timesteps/plt286093')
f_1600 = yt.load('/home/ltt/Output/Ign_Outflow/noout_methane01000')
f_1800 = yt.load('/home/ltt/Output/Ign_Outflow/noout_methane01500')

# Get the information needed for the grid
max_level = f_1600.index.max_level 

lo =  np.array([0.0, 0.0, 0])
hi = np.array([15.0, 5.0, 0.1])

dxmin = f_1600.index.get_smallest_dx()
dxmax = dxmin*2.0*2.0 # if amr level 2
npts=np.floor((hi-lo)/dxmin)
# lo=np.array([0.0,0.0,0.00])
# hi=np.array([7.2,1.8,0.45])
# lo=np.array([1.25,0.0,0.15])
# hi=np.array([1.65,0.4,0.35])
# lo=np.array([0.75,0.0,0.15])
# hi=np.array([1.35,0.4,0.35])
# dxmin = ds_286093.index.get_smallest_dx()
# dxmax = dxmin*2.0*2.0
# dx_np = np.array(dxmin)
# npts=np.floor((hi-lo)/dxmin)
# # npts=np.floor((hi-lo)/dxmax)
npts_np = np.array(npts)
# fields_load=["vfrac","density","pressure", "Temp", "x_velocity", "y_velocity", "z_velocity"]
# fields_load=["density","pressure", "Temp", "x_velocity", "y_velocity", "z_velocity",
#     "Y(N2)", "Y(H2)", "Y(H)", "Y(O)", "Y(O2)", "Y(H2O)", "Y(HO2)", "Y(H2O2)"]
fields_load=["density","pressure", "Temp", "x_velocity", "y_velocity", "z_velocity","Y(H2O)"]

# Get the values on the grid at level 0 (finest level)
first = f_1600.covering_grid(level=max_level, left_edge=lo, dims=npts, fields=fields_load)
second = f_1800.covering_grid(level=max_level, left_edge=lo, dims=npts, fields=fields_load)
# ad_286092 = ds_286092.covering_grid(level=0, left_edge=lo, dims=npts, fields=fields_load)
# ad_286093 = ds_286093.covering_grid(level=0, left_edge=lo, dims=npts, fields=fields_load)
pres_1600    = np.array(first["pressure"])
yh20         = np.array(first["Y(H2O)"])
pres_1800    = np.array(second["pressure"])


x_arr = np.linspace(lo[0], hi[0], int(npts_np[0]))
y_arr = np.linspace(lo[1], hi[1], int(npts_np[1]))

# Build time derivatives
# time_plt286092 = 0.0020000052369055794
# time_plt286093 = 0.0020000080533389458
time_plt1600 = 1.8331263994279536e-05
time_plt1800 = 2.8925491834175558e-05
dt = time_plt1800 - time_plt1600

check = pres_1600[:,:,0]
check2 = pres_1800[:,:,0]
max_index1 = np.unravel_index(np.argmax(check), check.shape)
max_index2 = np.unravel_index(np.argmax(check2), check2.shape)

speed = ((max_index2[0] - max_index1[0]) * dxmin) / (time_plt1800 - time_plt1600) *(10**-2)

print(speed, " m/s")

pressure = pres_1600[:,:,0]
pressure_2 = pres_1800[:,:,0]

contour = plt.contourf(x_arr,y_arr, pressure.transpose(), levels = 100, cmap="jet")

plt.colorbar(contour)

plt.savefig(f"pressure_1600.png")
plt.close()

contour = plt.contourf(x_arr,y_arr, pressure_2.transpose(), levels = 100, cmap="jet")

plt.colorbar(contour)

plt.savefig(f"pressure_1800.png")
plt.close()

#plt.plot(yh20[35,:,0].transpose(),pres_1600[35,:,0].transpose())
plt.plot(y_arr, pres_1600[35,:,0].transpose())
plt.savefig(f"yh20vpres_1600.png")
plt.close()




# ax.set_ylabel('Enthalpy on grid points')
# ax.legend()
# plt.savefig(f"full_case_enth.png")
# plt.close()