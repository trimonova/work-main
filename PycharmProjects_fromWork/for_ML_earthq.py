import numpy as np

import scipy.io as sio
import pandas as pd
import h5py
import os

total_list = []

main_folder = 'C:\disk_D\землетрясен. штуки\\2019-07-06\\2019-07-06\\'
files_list = os.listdir(main_folder)

for file in files_list:
    #file_number = line_list[1]
    #os.mkdir('C:\disk_D\землетрясен. штуки\\2019-07-06\\N_'+str(file_number))
    #f = scipy.io.loadmat('C:\disk_D\землетрясен. штуки\\2019-07-06\\2019-07-06\\N_'+str(file_number)+'.mat')
    #print(main_folder + file)
    f = h5py.File(main_folder + file, 'r')
    keys_list = list(f.keys())
    #print(keys_list)
    full_data = list((f.get('y')))
    #print(np.shape(full_data))
    if file != 'N_32.mat':
        for minute_number in range(1, 31):
            for setup_number in range(3):
                data_line = full_data[setup_number][(minute_number - 1) * 60 * 10000:minute_number * 60 * 10000]
                total_list.append({'file_name': file, 'minute_number': minute_number,
                                   'setup_number': setup_number, 'data_record': [data_line], 'time_in': [[]],
                                   'type': 0})
    else:
        for minute_number in range(1, 15):
            for setup_number in range(3):
                data_line = full_data[setup_number][(minute_number - 1) * 60 * 10000:minute_number * 60 * 10000]
                total_list.append({'file_name': file, 'minute_number': minute_number,
                                   'setup_number': setup_number, 'data_record': [data_line], 'time_in': [[]],
                                   'type': 0})
                # sio.savemat('C:\disk_D\землетрясен. штуки\\2019-07-06\\y.mat', {'y': data_line})

f = open('C:\disk_D\землетрясен. штуки\T3_int_t1t2t3t4.dat')
for line in f:
    line_list = line.split()
    if int(line_list[0]) > 10 and int(line_list[2] != 0):
        file_name = 'N_'+ line_list[1] + '.mat'

        for setup_number in range(3):
            time_in = float(line_list[setup_number + 2])
            time_in_min = int(float(line_list[setup_number + 2])/60)
            for elem in total_list:
                if (elem['file_name'] == file_name) and (elem['minute_number'] == time_in_min) and (elem['setup_number'] == setup_number):
                    elem['time_in'].append(round((time_in-time_in_min*60), 4))
                    elem['type'] = 2

#print(total_list)
for elem in total_list:
    if (elem['file_name'] == 'N_13.mat') and (elem['minute_number'] == 7) and (
            elem['setup_number'] == 0):
        print(elem['time_in'])

for elem in total_list:

    if ((elem['file_name'] == 'N_02.mat') or (elem['file_name'] == 'N_27.mat') or (elem['file_name'] == 'N_28.mat') or (elem['file_name'] == 'N_29.mat') or (elem['file_name'] == 'N_30.mat') or (elem['file_name'] == 'N_31.mat') or (elem['file_name'] == 'N_32.mat')):
        print('ok')
        elem['type'] = 1
print(total_list)

total_pd = pd.DataFrame()
for elem in total_list:
    elem_pd = pd.DataFrame(elem)
    total_pd.append(elem_pd)
print(total_pd)











