from src.CCMatrix import CCMatrix
import numpy as np


class Model:
    def __init__(self, cc, infectious_rate, incubation_rate, recovery_rate):
        # rates for model
        self.infectious_rate = infectious_rate
        self.incubation_rate = incubation_rate
        self.recovery_rate = recovery_rate
        self.contact_matrix = cc.cc_matrix * infectious_rate

        # initial data for model
        self.susceptible = cc.belgium_count.astype('float64')
        self.exposed = np.zeros(86)
        self.infected = np.zeros(86)
        self.infected[38] = 1
        self.recovered = np.zeros(86)

        # data (for graphing)
        self.infected_data = []
        self.recovered_data = []

    def infect(self):
        susceptible_transpose = np.transpose(np.asmatrix(self.susceptible))
        susceptible_infected = np.matmul(susceptible_transpose, np.asmatrix(self.infected))
        contacts = np.multiply(self.contact_matrix, susceptible_infected)
        contacts = np.asarray(contacts)

        for column in contacts.transpose():
            column = np.minimum(self.susceptible, column)
            self.infected += column
            self.susceptible -= column
        self.susceptible = np.maximum(self.susceptible, np.zeros(86))

    def recover(self, infected):
        recoveries = infected * self.recovery_rate
        self.recovered += recoveries
        self.infected -= recoveries
        self.infected = np.maximum(self.infected, np.zeros(86))

    def run(self, days):
        for i in range(days):
            infected = self.infected
            self.infect()
            self.recover(infected)

            self.infected_data.append(self.infected.sum())
            self.recovered_data.append(self.recovered.sum())
            # print(self.susceptible.sum() + infected.sum() + self.recovered.sum())
            # print(i)
        # print('---')
        # print('Number of susceptible: ', self.susceptible.sum())
        # print('Number of infected: ', self.infected.sum())
        # print('Number of recovered: ', self.recovered.sum())

        # print(self.susceptible.sum() + self.infected.sum() + self.recovered.sum())
