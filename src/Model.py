from src.CCMatrix import CCMatrix
import numpy as np


# TODO: aparte infected box voor gehospitaliseerden
# TODO: aparte susceptible doos voor mensen die in zorgsector werken (hogere contact rate dan de rest) ?
# TODO: lagere recovery rate voor gehospitaliseerden, ook hogere death rate
# TODO: eigenlijk gewoon model maken zoals op tekening

class Model:
    def __init__(self, cc, infectious_rate, incubation_rate, recovery_rate):
        # rates for model
        self.infectious_rate = infectious_rate
        self.incubation_rate = incubation_rate
        self.recovery_rate = recovery_rate
        self.contact_matrix = cc.cc_matrix * infectious_rate

        # initial data for model
        # SEIR
        self.susceptible = cc.belgium_count.astype('float64')
        self.exposed = np.zeros(86)
        self.infected = np.zeros(86)
        self.infected[38] = 1
        self.recovered = np.zeros(86)

        # Extra Compartments
        self.hospital = np.zeros(86)
        self.hospital_ic = np.zeros(86)

        # data (for graphing)
        self.infected_data = []
        self.recovered_data = []

    def infect(self):
        """
        From normal population to infected compartment
        :return:
        """
        susceptible_transpose = np.transpose(np.asmatrix(self.susceptible))
        # Normal Contact matrix for infected people
        susceptible_infected = np.matmul(susceptible_transpose, np.asmatrix(self.infected))
        # Lower Contact matrix for hospitalized people
        # TODO
        # Add both to get total infected
        # TODO

        contacts = np.multiply(self.contact_matrix, susceptible_infected)
        contacts = np.asarray(contacts)

        for column in contacts.transpose():
            column = np.minimum(self.susceptible, column)
            self.infected += column
            self.susceptible -= column
        self.susceptible = np.maximum(self.susceptible, np.zeros(86))

    def recover(self, infected):
        """
        From infected, hospital, ic to recovered compartment
        :param infected: infected people the day before
        :return:
        """
        recoveries = infected * self.recovery_rate
        self.recovered += recoveries
        self.infected -= recoveries
        self.infected = np.maximum(self.infected, np.zeros(86))

    def go_to_hospital(self, infected):
        """
        From infected compartment to hospital compartment
        :param infected: infected people the day before
        :return:
        """
        ...

    def go_to_ic(self, hospitalized):
        """
        From hospital compartment to IC compartment
        :param hospitalized: hospitalized people the day before
        :return:
        """
        ...

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
