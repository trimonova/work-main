import os
import pandas as pd
df = pd.read_excel('C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\exp 2\\xp_0_yp_0.xlsx')
print(df)
dir_name = 'C:\disk_D\Repos\PycharmProjects_fromWork\PycharmProjects\complicated_flow_model\\replacementModelAndWell\\first_exp\exp 2'
directory = os.fsencode(dir_name)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    print(filename)
    df = pd.read_excel(dir_name+'\\'+filename)
    print(df.iloc[:,1])
