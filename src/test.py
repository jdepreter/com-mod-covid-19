import pandas as pd
import numpy as np
from os.path import dirname, abspath
import os

data_folder = os.path.join(dirname(dirname(abspath(__file__))), 'data')
# Contact Matrix
matrix_csv = os.path.join(data_folder, 'COVID19BE_HOSP.csv')
df = pd.read_csv(matrix_csv, delimiter=',', )
df.groupby(["DATE"]).sum()
result = np.where(df['TOTAL_IN'] < 1, df['TOTAL_IN'], df['TOTAL_IN_ICU']/df['TOTAL_IN']).sum()/len(df)
print(result)
print("exit")
