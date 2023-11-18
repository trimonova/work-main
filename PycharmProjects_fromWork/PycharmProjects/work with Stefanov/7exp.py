import openpyxl
import pandas as pd
from matplotlib import pyplot as plt
from plot_graphs import plot_graphs
wb = openpyxl.load_workbook(filename = 'C:\\disk_D\work with Stefanov\\experimental data mathing to 0 (version 1).xlsx', data_only=True)
sheet = wb['UCS 7 эксп']

sigA = [v[0].value for v in sheet['E9:E578']]
epsA = [v[0].value for v in sheet['R9:R578']]
epsR = [v[0].value for v in sheet['S9:S578']]
epsV = [v[0].value for v in sheet['T9:T578']]

plt.plot(epsA, sigA)
plt.plot(epsR, sigA)
plt.plot(epsV, sigA)
plt.show()

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\20exp\\3D_Rock_V10_20exp_5\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


# plt.plot(epsA_calc, sigA_calc, 'r')
# plt.plot(epsR_calc, sigA_calc)
# plt.plot(epsV_calc, sigA_calc)
# plt.plot(epsA, sigA, color = 'r', linestyle = '--')
# plt.plot(epsR, sigA)
# plt.plot(epsV, sigA)
#
# plt.show()

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_6\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
                color_1='r', color_2='g', color_3='b', dot_style='--', label='6')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_8\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
                color_1='r', color_2='g', color_3='b', dot_style="-.", label='8')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\2expUCS\\3D_Rock_V10_2expUCS_15\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
                color_1='r', color_2='g', color_3='b', dot_style=":", label='15')

plt.legend(loc='lower left', prop={'size': 10})
plt.title('7 exp: 0.1 MPa')
plt.xlabel('Деформации')
plt.ylabel('Осевое напряжение, МПа')
plt.xlim((-0.01, 0.015))
plt.show()
# plt.plot(epsA_calc, sigA_calc, 'r')
# plt.plot(epsR_calc, sigA_calc)
# plt.plot(epsV_calc, sigA_calc)
# plt.plot(epsA, sigA, color = 'r', linestyle = '--')
# plt.plot(epsR, sigA)
# plt.plot(epsV, sigA)
#
# plt.show()
