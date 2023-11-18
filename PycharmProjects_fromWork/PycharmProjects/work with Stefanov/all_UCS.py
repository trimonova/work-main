import openpyxl
import pandas as pd
from matplotlib import pyplot as plt
from plot_graphs import plot_graphs

wb = openpyxl.load_workbook(filename = 'C:\\disk_D\work with Stefanov\\experimental data mathing to 0 (version 1).xlsx', data_only=True)
# sheet_2 = wb['UCS 2 эксп']
#
# sigA_2 = [v[0].value for v in sheet_2['E9:E482']]
# epsA_2 = [v[0].value for v in sheet_2['R9:R482']]
# epsR_2 = [v[0].value for v in sheet_2['S9:S482']]
# epsV_2 = [v[0].value for v in sheet_2['T9:T482']]
#
# plt.plot(epsA_2, sigA_2, color='r', linestyle='-', label='UCS 2')
# plt.plot(epsR_2, sigA_2, color='r', linestyle=':')
# plt.plot(epsV_2, sigA_2, color='r', linestyle='-.')
#
# sheet_6 = wb['UCS 6 эксп']
#
# sigA_6 = [v[0].value for v in sheet_6['E9:E360']]
# epsA_6 = [v[0].value for v in sheet_6['R9:R360']]
# epsR_6 = [v[0].value for v in sheet_6['S9:S360']]
# epsV_6 = [v[0].value for v in sheet_6['T9:T360']]
#
# plt.plot(epsA_6, sigA_6, color='g', linestyle='-', label='UCS 6')
# plt.plot(epsR_6, sigA_6, color='g', linestyle=':')
# plt.plot(epsV_6, sigA_6, color='g', linestyle='-.')

sheet_7 = wb['UCS 7 эксп']

sigA_7 = [v[0].value for v in sheet_7['E9:E578']]
epsA_7 = [v[0].value for v in sheet_7['R9:R578']]
epsR_7 = [v[0].value for v in sheet_7['S9:S578']]
epsV_7 = [v[0].value for v in sheet_7['T9:T578']]

plt.plot(epsA_7, sigA_7, color='b', linestyle='-', label='UCS 7')
plt.plot(epsR_7, sigA_7, color='b', linestyle=':')
plt.plot(epsV_7, sigA_7, color='b', linestyle='-.')

# sheet_8 = wb['UCS 8 эксп']
#
# sigA_8 = [v[0].value for v in sheet_8['E9:E629']]
# epsA_8 = [v[0].value for v in sheet_8['R9:R629']]
# epsR_8 = [v[0].value for v in sheet_8['S9:S629']]
# epsV_8 = [v[0].value for v in sheet_8['T9:T629']]
#
# plt.plot(epsA_8, sigA_8, color='y', linestyle='-', label='UCS 8')
# plt.plot(epsR_8, sigA_8, color='y', linestyle=':')
# plt.plot(epsV_8, sigA_8, color='y', linestyle='-.')

