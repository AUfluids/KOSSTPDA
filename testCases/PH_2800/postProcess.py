import numpy as np
import pandas as pd
import os
import Ofpp
print(os.getcwd())
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib as mlp
from matplotlib.lines import Line2D

from scipy.interpolate import griddata


plt.rcParams.update({
    "text.usetex": False,
    # 'text.latex.preamble': r"\usepackage{lmodern}",
    # 'font.family': 'lmodern',
    'font.size': 7,
    'axes.linewidth': 0.1,
    "grid.color": "#cccccc",
    "axes.grid": False,
    "grid.alpha": 0.8,
    "grid.linewidth": 0.5,
    "grid.linestyle": '-',
    "axes.grid.which": 'both',
    "axes.spines.right": False,
    "axes.spines.top": False,
    'axes.axisbelow': True})

# plt.rcParams.update({'font.family': 'Times New Roman',
#                      'mathtext.fontset': 'cm',
#                      'font.weight': 'normal',
#                      'font.size': 12,
#                     })


def VolumetricAverage(P_HF,Volume):
    return np.sum((np.multiply(np.squeeze(P_HF),Volume)))/np.sum(Volume)



C_name = './!RefData/C'
V_name = './!RefData/V'
C_ph = Ofpp.parse_internal_field(C_name)
V_ph = Ofpp.parse_internal_field(V_name)
C_x, C_y, C_z = np.split(C_ph,3,axis=1)

C_boundary = Ofpp.parse_boundary_field(C_name)
bottomWall = C_boundary[b'bottomWall'][b'value'] 
gap = 0.03

def plotFrame(bottomWall, gap):
    plt.plot(bottomWall[:,0],bottomWall[:,1],'k')
    plt.plot([0+gap,0+gap],[1,3.035],'k')
    plt.plot([9-gap,9-gap],[1,3.035],'k')
    plt.plot([0+gap,9-gap],[3.035,3.035],'k')
    plt.xticks([2,4,6,8], ['2','4','6','8'])
    plt.xlabel(r"$x/H$")
    plt.yticks([0,1,2,3], ['0','1','2','3'])
    plt.ylabel(r"$y/H$")
    plt.xlim([0,9])
    plt.ylim([0,3.035])


C_x = np.squeeze(C_x)
C_y = np.squeeze(C_y)



# High Fidelity data


fname = './!RefData/U_HF'
U_internal = Ofpp.parse_internal_field(fname)
U_x_LES, U_y_LES, U_z_LES = np.split(U_internal,3,axis=1)
U_mag_LES = np.sqrt(U_x_LES**2 + U_y_LES**2 + U_z_LES**2)



Cf_HF = np.genfromtxt('./!RefData/CF_HF.txt')
print(Cf_HF.shape)


# U_balk = VolumetricAverage(U_x_LES,V_ph)
U_balk = 1
print('Ubalk is ', U_balk)

bottomWS_LES = griddata(Cf_HF[:,0], Cf_HF[:,1], bottomWall[:,0])*((0.5*U_balk*U_balk))


# fname = './!RefData/wallShearStress_HF'
# wallShear = Ofpp.parse_boundary_field(fname)
# bottomWallShear = wallShear[b'bottomWall'][b'value'] 
# bottomST_x, bottomST_y, bottomST_z = np.split(bottomWallShear,3,axis=1)
# bottomWS_LES = np.sqrt(np.square(bottomST_x) + np.square(bottomST_y) + np.square(bottomST_z))*(-bottomST_x/np.abs(bottomST_x))


# Baseline data

fname = './!RefData/U_KOSST'
U_internal = Ofpp.parse_internal_field(fname)
U_x_KOSST, U_y_KOSST, U_z_KOSST = np.split(U_internal,3,axis=1)
U_mag_KOSST = np.sqrt(U_x_KOSST**2 + U_y_KOSST**2 + U_z_KOSST**2)


fname = './!RefData/wallShearStress_KOSST'
wallShear = Ofpp.parse_boundary_field(fname)
bottomWallShear = wallShear[b'bottomWall'][b'value'] 
bottomST_x, bottomST_y, bottomST_z = np.split(bottomWallShear,3,axis=1)
bottomWS_KOSST = np.sqrt(np.square(bottomST_x) + np.square(bottomST_y) + np.square(bottomST_z))*(-bottomST_x/np.abs(bottomST_x))


# Reading this case

list_subfolders = [f.name for f in os.scandir(os.getcwd()) if f.is_dir()]
print(list_subfolders)

list_subfoders_num = []
for folder in list_subfolders:
    if folder != '0':
        if (folder.isnumeric()):
            list_subfoders_num.append(folder)
list_subfoders_num.sort()
res = [eval(i) for i in list_subfoders_num]
selectedFolder = str(np.max(res))
print(selectedFolder)



fname = os.getcwd()+ '/' + selectedFolder + '/U'
U_internal = Ofpp.parse_internal_field(fname)
U_x_, U_y_, U_z_ = np.split(U_internal,3,axis=1)
U_mag_ = np.sqrt(U_x_**2 + U_y_**2 + U_z_**2)


fname = os.getcwd()+ '/' + selectedFolder + '/wallShearStress'
wallShear = Ofpp.parse_boundary_field(fname)
bottomWallShear = wallShear[b'bottomWall'][b'value'] 
bottomST_x, bottomST_y, bottomST_z = np.split(bottomWallShear,3,axis=1)
bottomWS_ = np.sqrt(np.square(bottomST_x) + np.square(bottomST_y) + np.square(bottomST_z))*(-bottomST_x/np.abs(bottomST_x))


# Plotting the contourplot with separation


plt.clf()
U_x_sqz = np.squeeze(U_x_)
U_y_sqz = np.squeeze(U_y_)
U_mag_sqz = np.sqrt(U_x_sqz**2 + U_y_sqz**2)/U_balk
U_mag_sqz_SST = np.sqrt(np.squeeze(U_x_KOSST)**2 + np.squeeze(U_y_KOSST)**2)/U_balk
U_mag_sqz_LES = np.sqrt(np.squeeze(U_x_LES)**2 + np.squeeze(U_y_LES)**2)/U_balk

U_mag_err = (U_mag_sqz - U_mag_sqz_LES) / VolumetricAverage(np.abs((U_mag_sqz_SST - U_mag_sqz_LES)), V_ph)

print('Minimum error: ' + str(U_mag_err.min()))
print('Maximum error: ' + str(U_mag_err.max()))

cmap = 'bwr'
vmin = -0.5
vmax = 1.5

cmapErr = 'RdBu_r'
vminErr = -1
vmaxErr = 1

triang = tri.Triangulation(C_x,C_y)
max_radius = 0.3
triangles = triang.triangles
# Mask off unwanted triangles.
xtri = C_x[triangles] - np.roll(C_x[triangles], 1, axis=1)
ytri = C_y[triangles] - np.roll(C_y[triangles], 1, axis=1)
maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
triang.set_mask(maxi > max_radius)


# plt.tricontourf(triang,
#                 U_x_sqz/U_balk,
#                 50,
#                 cmap=cmap,
#                 vmin=vmin, 
#                 vmax=vmax
#                 )
plt.tricontourf(triang,
                U_mag_err,
                50,
                cmap=cmapErr,
                vmin=vminErr, 
                vmax=vmaxErr
                )


xi = np.linspace(C_x.min(), C_x.max(), 50)
yi = np.linspace(C_y.min(), C_y.max(), 50)
X, Y = np.meshgrid(xi, yi)
U = interpolate.griddata((C_x, C_y), U_x_sqz, (X, Y), method='cubic')
V = interpolate.griddata((C_x, C_y), U_y_sqz, (X, Y), method='cubic')

plt.streamplot(X, Y, U, V, linewidth=0.5, color= 'k', density=1, arrowsize=0.8, arrowstyle='-|>')

plotFrame(bottomWall, gap)

NormErrorUmag = VolumetricAverage(np.abs(U_mag_ - U_mag_LES), V_ph)/VolumetricAverage(np.abs(U_mag_KOSST - U_mag_LES), V_ph)

NormErrorWSS = np.mean(np.abs(np.squeeze(bottomWS_) - bottomWS_LES))/np.mean(np.abs(np.squeeze(bottomWS_KOSST) - bottomWS_LES))

NormErrorAvg = (NormErrorUmag + NormErrorWSS) / 2

# print(bottomWS_.shape)
# print(bottomWS_KOSST.shape)
# print(bottomWS_LES.shape)

# print(np.abs(bottomWS_ - bottomWS_LES))
# print(np.abs(bottomWS_KOSST - bottomWS_LES))


if selectedFolder == '3500':
    if NormErrorUmag <= 2:
        NormErrorUmag = 2
    if NormErrorWSS <= 2:
        NormErrorWSS = 2
    NormErrorAvg = (NormErrorUmag + NormErrorWSS) / 2

fileName = os.getcwd() +'/objFunResults.txt'
with open(fileName, 'w') as f:
        f.write('j1 = ')
        f.write('%.8f' %(NormErrorUmag))
        f.write('\nj2 = ')
        f.write('%.8f' %(NormErrorWSS))
        f.write('\nJ = ')
        f.write('%.8f' %(NormErrorAvg))


plt.title('$J = ' + str(np.round(NormErrorAvg, 3)) + ', j_{1} = ' + str(np.round(NormErrorUmag, 3)) + ', j_{2} = ' + str(np.round(NormErrorWSS, 3)) + '$',
    loc='left')

ax = plt.gca() #you first need to get the axis handle
ax.set_aspect(1) #sets the height to width ratio to 1.5.
norm = mlp.colors.Normalize(vmin=vmin, vmax=vmax)
normErr = mlp.colors.Normalize(vmin=vminErr, vmax=vmaxErr)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right",
                            size="2%",
                            pad=0.1)
# cbar = mlp.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,label = r"$(\langle u \rangle^2 + \langle v \rangle^2)^{0.5}/u_{b}$")
# cbar.set_ticks([vmin, 0.8, vmax]) #[0, 0.01, 0.03, 0.05] [0, 0.5, 1, 1.5] [0, 5, 10, 15] [0, 0.3, 0.7, 1]
# cbar.set_ticklabels([vmin, 0.8, vmax])

cbar = mlp.colorbar.ColorbarBase(cax, cmap=cmapErr, norm=normErr,label = r"$\frac{\| \langle u \rangle \| - \| \langle u_{LES} \rangle \|}{\| \langle u_{SST} \rangle \| - \| \langle u_{LES} \rangle \|}$")
cbar.set_ticks([vminErr, -0.5, 0, 0.5, vmaxErr]) #[0, 0.01, 0.03, 0.05] [0, 0.5, 1, 1.5] [0, 5, 10, 15] [0, 0.3, 0.7, 1]
cbar.set_ticklabels([vminErr, -0.5, 0, 0.5, vmaxErr])


nameTosave = 'plt_Stream.png'
plt.savefig(nameTosave, dpi=300, bbox_inches='tight')




#  Plotting the friction coefficient

plt.clf()

plt.rcParams.update({
    "text.usetex": False,
    # 'text.latex.preamble': r"\usepackage{lmodern}",
    # 'font.family': 'lmodern',
    'font.size': 7,
    'axes.linewidth': 2,
    "grid.color": "#cccccc",
    "axes.grid": True,
    "grid.alpha": 0.8,
    "grid.linewidth": 0.5,
    "grid.linestyle": '-',
    "axes.grid.which": 'both',
    "axes.spines.right": False,
    "axes.spines.top": False,
    'axes.axisbelow': True})

plt.plot(bottomWall[:,0],bottomWS_LES/(0.5*U_balk*U_balk),color = 'k', marker = 'o', markersize = 3, linestyle = 'none' , label = 'HF')
plt.plot(bottomWall[:,0],bottomWS_KOSST/(0.5*U_balk*U_balk),'-g', label = 'KOSST')
plt.plot(bottomWall[:,0],bottomWS_/(0.5*U_balk*U_balk),'--r', label = 'This case')

# plt.title('$J = ' + str(np.round(NormErrorAvg, 3)) + '$J_{mag} = ' + str(np.round(NormErrorUmag, 3)) + ', J_{tau} = ' + str(np.round(NormErrorWSS, 3)) + ', J_{tke} = ' + str(np.round(NormErrorTKE, 3)) + '$',loc='left')
plt.title('$J = ' + str(np.round(NormErrorAvg, 3)) + ', j_{1} = ' + str(np.round(NormErrorUmag, 3)) + ', j_{2} = ' + str(np.round(NormErrorWSS, 3)) + '$',
    loc='left')

# plt.xticks([2,4,6,8], ['2','4','6','8'])
plt.xlabel(r"$x/H$")
# plt.yticks([0,1,2,3], ['0','1','2','3'])
plt.ylabel(r"$C_f$")
plt.legend(loc = 'best')
nameTosave = 'plt_Cf.png'
plt.savefig(nameTosave, dpi=300, bbox_inches='tight')



