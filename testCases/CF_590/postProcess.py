import numpy as np
import os
import Ofpp
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams.update({
    "text.usetex": False,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
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

def VolumetricAverage(P_HF,Volume):
    return np.sum((np.multiply(np.squeeze(P_HF),Volume))) / np.sum(Volume)


cm = 1 / 2.54  # inches to cm
textWidth = 9  # cm
figSize = (textWidth * cm, textWidth * cm * 3 / 4)

figSizeMedium = (14 * cm, 14 * cm * 3 / 4)
figSizeFull = (19 * cm, 19 * cm)

ReTau = 590

C_name = './!RefData/C'
V_name = './!RefData/V'
C_y = Ofpp.parse_internal_field(C_name)
V_CF = Ofpp.parse_internal_field(V_name)

# print(C_y)

# C_x, C_y, C_z = np.split(C_CF,3,axis=1)

HF_Cf = np.genfromtxt('./!RefData/chan590.means.txt')  # PH_Martin2018.csv  TwallBruer
HF_Cf_Rij = np.genfromtxt('./!RefData/chan590.reystress.txt')  # PH_Martin2018.csv  TwallBruer
k_HF = (HF_Cf_Rij[:, 2] + HF_Cf_Rij[:, 3] + HF_Cf_Rij[:, 4]) / 2

fname = './!RefData/U_KOSST'
U_internal_ = Ofpp.parse_internal_field(fname)
U_x_SST, U_y_SST, U_z_SST = np.split(U_internal_, 3, axis=1)

fname = './!RefData/k_KOSST'
k_SST = Ofpp.parse_internal_field(fname)


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


fname = './!RefData/C'
y_values = Ofpp.parse_internal_field(fname)
y_plus_KOSST = y_values * ReTau

fname = './!RefData/C'
y_values = Ofpp.parse_internal_field(fname)
y_plus = y_values * ReTau

fname = selectedFolder + '/U'
U_internal_ = Ofpp.parse_internal_field(fname)
U_x_, U_y_, U_z_ = np.split(U_internal_,3,axis=1)

fname = selectedFolder + '/k'
k_ = Ofpp.parse_internal_field(fname)

U_x_err = U_x_ - U_x_SST

# print(V_CF)

U_x_err_vol = VolumetricAverage(np.abs(U_x_ - U_x_SST), V_CF)
U_x_vol_KOSST = VolumetricAverage(U_x_SST, V_CF)

U_x_err_abs = U_x_err_vol / U_x_vol_KOSST  # VolumetricAverage(np.abs(U_x_ - U_x_SST), V_CF) / VolumetricAverage(U_x_SST, V_CF)

k_err_vol = VolumetricAverage(np.abs(k_ - k_SST), V_CF)
k_vol_KOSST = VolumetricAverage(k_SST, V_CF)

k_err_abs = k_err_vol / k_vol_KOSST  # VolumetricAverage(np.abs(U_x_ - U_x_SST), V_CF) / VolumetricAverage(U_x_SST, V_CF)

print('Error of u:' + str(U_x_err_abs))
print('Error of TKE:' + str(k_err_abs))

fig, ax = plt.subplots(figsize=figSizeMedium, ncols=2, nrows=1, sharey=False)

ax[0].plot(HF_Cf[:, 1],HF_Cf[:, 2], linestyle='None', marker='o',
            markerfacecolor='None',
            markeredgecolor='k', 
            markersize=5,
            label='DNS ' + str(ReTau))
ax[0].plot(y_plus_KOSST, U_x_SST ,'-', c='tab:green', label = r"$k-\omega$ SST")
ax[0].plot(y_plus, U_x_,'--', c='tab:red', label = "This Case")

ax[0].set_xlabel(r"$y^+$")
ax[0].set_ylabel(r"$U^+$")
ax[0].legend(loc = 'best')
ax[0].set_xlim([0.1, ReTau])
ax[0].set_xscale("log")

ax[1].plot(HF_Cf[:,1],k_HF, linestyle='None', marker='o',
            markerfacecolor='None',
            markeredgecolor='k', 
            markersize=5,
            label='DNS ' + str(ReTau))
ax[1].plot(y_plus_KOSST, k_SST ,'-', c='tab:green', label = r"$k-\omega$ SST")
ax[1].plot(y_plus, k_,'--', c='tab:red', label = "This Case")

ax[1].set_xlabel(r"$y^+$")
ax[1].set_ylabel(r"$k$")
ax[1].set_xlim([0, ReTau])

nameTosave = 'Ukprofile'  + '.png'
fig.savefig(nameTosave, dpi=300, transparent=False)


# if float(selectedFolder) < 20000: 
j1 = U_x_err_abs
j2 = k_err_abs
J = (j1 + j2) / 2
# else:
#     j1 = 1
#     j2 = 1
#     J = 1

fileName = os.getcwd() +'/objFunResults.txt'
with open(fileName, 'w') as f:
        f.write('j1 = ')
        f.write('%.8f' %(j1))
        f.write('\nj2 = ')
        f.write('%.8f' %(j2))
        f.write('\nJ = ')
        f.write('%.8f' %(J))
