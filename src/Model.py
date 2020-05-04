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
        self.recovery_rate_hospital = recovery_rate / 2.0
        self.recovery_rate_ic = recovery_rate / 3.0
        self.contact_matrix = cc.cc_matrix * infectious_rate

        # extra rates
        self.hospital_rate = np.array([(85-i)/10 for i in range(86)])     # Temp value
        self.hospital_ic_rate = np.full(86, 0.2)  # Temp value
        self.death_rate = np.full(86, 0.26)        # Temp value

        # initial data for model
        # SEIR
        self.susceptible = cc.belgium_count.astype('float64')
        self.exposed = np.zeros(86)
        self.infected = np.zeros(86)
        self.infected[38] = 1
        self.recovered = np.zeros(86)

        # Extra Compartments
        self.susceptible_hospital_staff = np.append(np.zeros(21), np.append(np.full(45, 3600), np.zeros(20)))
        self.susceptible -= self.susceptible_hospital_staff
        self.hospital = np.zeros(86)
        self.hospital_ic = np.zeros(86)
        self.dead = np.zeros(86)

        # data (for graphing)
        self.infected_data = []
        self.recovered_data = []
        self.hospital_data = []
        self.ic_data = []
        self.dead_data = []

    def infect(self, susceptible, infected, factor):
        """
        From normal population to infected compartment
        :return:
        """
        susceptible_transpose = np.transpose(np.asmatrix(susceptible))
        # Normal Contact matrix for infected people
        susceptible_infected = np.matmul(susceptible_transpose, np.asmatrix(infected))

        contacts = np.multiply(self.contact_matrix, susceptible_infected) * factor
        contacts = np.asarray(contacts)

        for column in contacts.transpose():
            column = np.minimum(self.susceptible, column)
            self.infected += column
            self.susceptible -= column
        self.susceptible = np.maximum(self.susceptible, np.zeros(86))


    def recover(self, infected, hospitalized, ic):
        """
        From infected, hospital, ic to recovered compartment
        :param infected: infected people the day before
        :param hospitalized: hospitalized people the day before
        :param ic: people in ic the day before
        :return:
        """
        recoveries = infected * self.recovery_rate
        self.recovered += recoveries
        self.infected -= recoveries
        self.infected = np.maximum(self.infected, np.zeros(86))

        recoveries = hospitalized * self.recovery_rate_hospital
        self.recovered += recoveries
        self.hospital -= recoveries
        self.hospital = np.maximum(self.hospital, np.zeros(86))

        recoveries = ic * self.recovery_rate_ic
        self.recovered += recoveries
        self.hospital_ic -= recoveries
        self.hospital_ic = np.maximum(self.hospital_ic, np.zeros(86))

    def go_to_hospital(self, infected):
        """
        From infected compartment to hospital compartment
        :param infected: infected people the day before
        :return:
        """
        hospitalized = infected * self.hospital_rate
        self.hospital += hospitalized
        self.infected -= hospitalized
        self.infected = np.maximum(self.infected, np.zeros(86))

    def go_to_ic(self, hospitalized):
        """
        From hospital compartment to IC compartment
        :param hospitalized: hospitalized people the day before
        :return:
        """
        ic = hospitalized * self.hospital_ic_rate
        self.hospital_ic += ic
        self.hospital -= ic
        self.hospital = np.maximum(self.hospital, np.zeros(86))
        ...

    def die(self, infected, hospitalized, ic):
        """
        From infected, hospitalized, ic to Death
        :param infected: infected people (not in hospital / ic)
        :param hospitalized: people in hospital (not in ic)
        :param ic: people in ic
        :return:
        """
        dead = ic * self.death_rate
        self.dead += dead
        self.hospital_ic -= dead
        self.hospital_ic = np.maximum(self.hospital_ic, np.zeros(86))
        # TODO infected -> dead, hospital -> dead
        ...

    def run(self, days):
        for i in range(days):
            # Do not use new data for calculating transitions
            infected = self.infected
            hospital = self.hospital
            ic = self.hospital_ic

            # Transitions between compartments
            self.infect(self.susceptible, infected, 1)
            self.infect(self.susceptible_hospital_staff, infected, 5)
            self.recover(infected, hospital, ic)
            self.go_to_hospital(infected)
            self.go_to_ic(hospital)
            self.die(infected, hospital, ic)

            # Save data for graphing
            self.infected_data.append(self.infected.sum())
            self.recovered_data.append(self.recovered.sum())
            self.hospital_data.append(self.hospital.sum())
            self.ic_data.append(self.hospital_ic.sum())
            self.dead_data.append(self.dead.sum())

            # print(self.susceptible.sum() + infected.sum() + self.recovered.sum())
            # print(i)
        # print('---')
        # print('Number of susceptible: ', self.susceptible.sum())
        # print('Number of infected: ', self.infected.sum())
        # print('Number of recovered: ', self.recovered.sum())

        # print(self.susceptible.sum() + self.infected.sum() + self.recovered.sum())
