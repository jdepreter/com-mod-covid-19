from numpy import genfromtxt
import pandas as pd
from os.path import dirname, abspath
import os


class CCMatrix:
    def __init__(self, matrix_csv, ages_csv):
        data_folder = os.path.join(dirname(dirname(abspath(__file__))), 'data')

        # Contact Matrix
        matrix_csv = os.path.join(data_folder, matrix_csv)
        self.cc_matrix = pd.read_csv(matrix_csv, header=None, delimiter=';').to_numpy()

        # Load pop data
        belgium_italy_csv = os.path.join(data_folder, ages_csv)
        self.belgium_italy_df = pd.read_csv(belgium_italy_csv, delimiter=',', index_col=None)
        # Split data from Belgium & Italy
        belgian = self.belgium_italy_df["GEO"] == 'Belgium'
        italian = self.belgium_italy_df["GEO"] == 'Italy'
        # Only use ages <= 85
        age_ok = self.belgium_italy_df["AGE"] < 86
        self.belgium_count = self.belgium_italy_df[belgian & age_ok]
        self.belgium_count = self.belgium_count["Value"].to_numpy()     # Create numpy array
        self.italy_count = self.belgium_italy_df[italian & age_ok]
        self.italy_count = self.italy_count["Value"].to_numpy()         # Create numpy array

        self.belgium_distribution = self.belgium_count/self.belgium_count.sum()
        self.italy_distribution = self.italy_count/self.italy_count.sum()
        self.belgian_italy_count = self.belgium_distribution/self.italy_distribution*self.italy_count
        self.italic_belgium_count = self.italy_distribution/self.belgium_distribution*self.belgium_count

        print('cc matrix with size', self.cc_matrix.size, 'inserted')
        print('age array with size', self.belgium_count.size, 'inserted')


# c = CCMatrix('cc15.csv', 'eurostat_pop_age.csv')
