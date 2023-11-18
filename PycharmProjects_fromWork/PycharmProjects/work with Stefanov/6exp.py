import openpyxl
import pandas as pd
from matplotlib import pyplot as plt
wb = openpyxl.load_workbook(filename = 'C:\\disk_D\work with Stefanov\\experimental data mathing to 0 (version 1).xlsx', data_only=True)
sheet = wb['UCS 6 эксп']

sigA = [v[0].value for v in sheet['E9:E360']]
epsA = [v[0].value for v in sheet['R9:R360']]
epsR = [v[0].value for v in sheet['S9:S360']]
epsV = [v[0].value for v in sheet['T9:T360']]

plt.plot(epsA, sigA)
plt.plot(epsR, sigA)
plt.plot(epsV, sigA)
plt.show()

df = pd.read_table("C:\disk_D\work with Stefanov\\3D_Rock_V10\\20exp\\3D_Rock_V10_20exp_5\\Diagramm.dat", sep="\s+")
sigA_calc = df['-SigmaQ(MPa)']
epsA_calc = df['-Eps1*10^5']/10**5
epsR_calc = df['-EpsR*10^5']/10**5
epsV_calc = df['-dV*10^5']/10**5


plt.plot(epsA_calc, sigA_calc, 'r')
plt.plot(epsR_calc, sigA_calc)
plt.plot(epsV_calc, sigA_calc)
plt.plot(epsA, sigA, color = 'r', linestyle = '--')
plt.plot(epsR, sigA)
plt.plot(epsV, sigA)

plt.show()