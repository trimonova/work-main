import openpyxl
import pandas as pd
from matplotlib import pyplot as plt
from plot_graphs import plot_graphs
wb = openpyxl.load_workbook(filename = 'C:\\disk_D\work with Stefanov\\experimental data mathing to 0 (version 1).xlsx', data_only=True)
sheet = wb['Experiment 1 (19)']

sigA = [v[0].value for v in sheet['E9:E405']]
epsA = [v[0].value for v in sheet['O9:O405']]
epsR = [v[0].value for v in sheet['P9:P405']]
epsV = [v[0].value for v in sheet['Q9:Q405']]

plt.plot(epsA, sigA)
plt.plot(epsR, sigA)
plt.plot(epsV, sigA)
plt.show()

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_6\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
                color_1='r', color_2='g', color_3='b', dot_style='--', label='6')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_4\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                  color_1='r', color_2='g', color_3='b', dot_style=(0, (12, 10)), label='4')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_3\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style='-.', label='3')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_2\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style=(0, (6, 10)), label='2')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style=(0, (1, 10)), label='1')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_8\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
                color_1='r', color_2='g', color_3='b', dot_style="-.", label='8')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_7\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style=(0, (9, 10)), label='7')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_15\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5

# changes_file = open("C:\disk_D\work with Stefanov\\3D_Rock_V10\\20exp\\3D_Rock_V10_20exp\changes_for_py.txt", "r")
# changes = changes_file.readline()

# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style=':', label='15')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_9\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style='-.', label='9')

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\19exp\\3D_Rock_V10_19exp_11\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


# plot_graphs(epsA_calc, epsR_calc, epsV_calc, sigA_calc, epsA, epsR, epsV, sigA,
#                 color_1='r', color_2='g', color_3='b', dot_style='--', label='11')

plt.legend(loc='lower left', prop={'size': 10})
plt.title('19 exp: 0.7 MPa')
plt.xlabel('Деформации')
plt.ylabel('Осевое напряжение, МПа')
plt.xlim((-0.01, 0.015))
plt.show()