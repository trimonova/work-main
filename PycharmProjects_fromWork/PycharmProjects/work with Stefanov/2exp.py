import openpyxl
import pandas as pd
from matplotlib import pyplot as plt
from plot_graphs import plot_graphs
import all_UCS

# wb = openpyxl.load_workbook(filename = 'C:\\disk_D\work with Stefanov\\experimental data mathing to 0 (version 1).xlsx', data_only=True)
# sheet = wb['UCS 2 эксп']
#
# sigA = [v[0].value for v in sheet['E9:E482']]
# epsA = [v[0].value for v in sheet['R9:R482']]
# epsR = [v[0].value for v in sheet['S9:S482']]
# epsV = [v[0].value for v in sheet['T9:T482']]
#
# plt.plot(epsA, sigA)
# plt.plot(epsR, sigA)
# plt.plot(epsV, sigA)
# plt.show()

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_7\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

plt.plot(epsA_calc, sigA_calc, color='black', linewidth=3, linestyle='-', label='UCS calc 7')
plt.plot(epsR_calc, sigA_calc, color='black', linewidth=3, linestyle=':')
plt.plot(epsV_calc, sigA_calc, color='black', linewidth=3, linestyle='-.')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_6\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

plt.plot(epsA_calc, sigA_calc, color='grey', linewidth=3, linestyle='-', label='UCS calc 6')
plt.plot(epsR_calc, sigA_calc, color='grey', linewidth=3, linestyle=':')
plt.plot(epsV_calc, sigA_calc, color='grey', linewidth=3, linestyle='-.')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_5\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

plt.plot(epsA_calc, sigA_calc, color='orange', linewidth=3, linestyle='-', label='UCS calc 5')
plt.plot(epsR_calc, sigA_calc, color='orange', linewidth=3, linestyle=':')
plt.plot(epsV_calc, sigA_calc, color='orange', linewidth=3, linestyle='-.')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_8\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

plt.plot(epsA_calc, sigA_calc, color='violet', linewidth=3, linestyle='-', label='UCS calc 8')
plt.plot(epsR_calc, sigA_calc, color='violet', linewidth=3, linestyle=':')
plt.plot(epsV_calc, sigA_calc, color='violet', linewidth=3, linestyle='-.')

plt.legend()
plt.show()

# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style='--', label='5')
#
# plt.legend(loc='lower left')
# plt.title('2 exp: 0.1 MPa')
# plt.xlabel('Деформации')
# plt.ylabel('Осевое напряжение, МПа')
# plt.show()
