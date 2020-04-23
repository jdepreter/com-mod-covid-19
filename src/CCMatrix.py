from numpy import genfromtxt
import pandas as pd
from os.path import dirname, abspath
import os



class CCMatrix:
    def __init__(self, matrix_csv, ages_csv):
        data_folder = os.path.join(dirname(dirname(abspath(__file__))), 'data')
        matrix_csv = os.path.join(data_folder, matrix_csv)
        ages_csv = os.path.join(data_folder, ages_csv)
        self.cc_matrix = pd.read_csv(matrix_csv, header=None, delimiter=';').to_numpy()
        self.age_array = genfromtxt(ages_csv, delimiter=',')[:86]
        print('cc matrix with size', self.cc_matrix.size, 'inserted')
        print('age array with size', self.age_array.size, 'inserted')
