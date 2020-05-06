import pandas as pd
import numpy as np
from os.path import dirname, abspath
import os

data_folder = os.path.join(dirname(dirname(abspath(__file__))), 'com-mod-covid-19/data')
# Contact Matrix
matrix_csv = os.path.join(data_folder, 'COVID19BE_HOSP.csv')
covid_belgie = pd.read_csv(matrix_csv, delimiter=',', )
covid_belgie.groupby(["DATE"]).sum()

print("exit")
