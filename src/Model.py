from src.CCMatrix import CCMatrix
import numpy as np
from scipy import interpolate


class Model:
    def __init__(self, contact_matrix, susceptible, infectious_rate, measure_factor=1.0, measure_day=0,
                 scenario=None, first_patient_age=38, cap_ic=False):
        # infectious / incubation rate
        self.infectious_rate = infectious_rate
        self.incubation_rate = 1.0 / 3.0

        # chances for hospitalization / intensive care / death
        self.hospital_chance = 0.19
        self.ic_chance = 0.224
        self.hospital_death_chance = 0.04
        self.ic_death_chance = 0.26
        self.ic_no_bed_death_chance = 0.99

        # recovery rates
        self.recovery_rate = (1.0 - self.hospital_chance) / 6.0
        self.recovery_rate_hospital = (1.0 - self.ic_chance - self.hospital_death_chance) / 8.0
        self.recovery_rate_ic = (1.0 - self.ic_death_chance) / 10.0
        self.recovery_rate_ic_no_bed = (1.0 - self.ic_no_bed_death_chance) / 10.0

        # contact matrix (combined with infectious rate and possible measures)
        self.contact_matrix = contact_matrix * infectious_rate

        # rates for hospital / ic / death
        self.hospital_rate = self.hospital_chance / 6.0
        self.ic_rate = np.full(86, self.ic_chance / 8.0)
        self.death_rate = np.full(86, self.ic_death_chance / 10.0)
        self.hospital_death_rate = np.full(86, self.hospital_death_chance / 8.0)
        self.ic_no_bed_death_rate = np.full(86, self.ic_no_bed_death_chance / 10.0)

        # initial data for model
        self.susceptible = susceptible.astype('float64')
        self.exposed = np.zeros(86)
        self.infected = np.zeros(86)
        self.infected[first_patient_age] = 1
        self.recovered = np.zeros(86)

        # Extra Compartments
        self.hospital = np.zeros(86)
        self.hospital_total = np.zeros(86)
        self.hospital_ic = np.zeros(86)
        self.dead = np.zeros(86)

        # measures
        self.measure_factor = measure_factor
        self.measure_day = measure_day

        # icu capacity (3000, 75 for each age below 40)
        self.ic_capacity = np.array([75 if i < 40 else 0 for i in range(86)])
        self.cap_ic = cap_ic

        # data (for graphing)
        self.infected_data = np.empty(0)
        self.exposed_data = np.empty(0)
        self.recovered_data = np.empty(0)
        self.hospital_data = np.empty(0)
        self.hospital_total_data = np.empty(0)
        self.ic_data = np.empty(0)
        self.dead_data = np.empty(0)
        self.case_data = np.empty(0)
        self.scenario = scenario

    def infect(self, susceptible, infected, factor):
        """
        From normal population to infected compartment
        :return:
        """
        susceptible_transpose = np.transpose(np.asmatrix(susceptible))
        # Normal Contact matrix for infected people
        susceptible_infected = np.matmul(susceptible_transpose, np.asmatrix(infected))

        contacts = np.multiply(np.multiply(self.contact_matrix, susceptible_infected), factor)
        contacts = np.asarray(contacts)

        for column in contacts.transpose():
            column = np.minimum(self.susceptible, column)
            self.exposed += column
            self.susceptible -= column
        self.susceptible = np.maximum(self.susceptible, np.zeros(86))

    def exp_to_inf(self, exposed):
        exp_to_inf = self.incubation_rate * exposed
        self.infected += exp_to_inf
        self.exposed -= exp_to_inf

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

        # ic recovery rate is lower when
        recoveries = ic * self.recovery_rate_ic
        if self.cap_ic:
            ic_no_bed = np.maximum(ic - self.ic_capacity, np.zeros(86))
            # if icu capacity is exceeded, let people without a bed recover at a low rate
            if ic_no_bed.sum() > 0:
                recoveries = self.ic_capacity * self.recovery_rate_ic
                recoveries += ic_no_bed * self.recovery_rate_ic_no_bed
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
        self.hospital_total += hospitalized
        self.infected -= hospitalized
        self.infected = np.maximum(self.infected, np.zeros(86))

    def go_to_ic(self, hospitalized):
        """
        From hospital compartment to IC compartment
        :param hospitalized: hospitalized people the day before
        :return:
        """
        ic = hospitalized * self.ic_rate
        self.hospital_ic += ic
        self.hospital -= ic
        self.hospital = np.maximum(self.hospital, np.zeros(86))
        ...

    def die(self, hospitalized, ic):
        """
        From hospitalized and ic to Death
        :param infected: infected people (not in hospital / ic)
        :param hospitalized: people in hospital (not in ic)
        :param ic: people in ic
        :return:
        """
        # icu dead
        dead = ic * self.death_rate
        if self.cap_ic:
            ic_no_bed = np.maximum(ic - self.ic_capacity, np.zeros(86))
            # if icu capacity is exceeded, let people without a bed die with at a high rate
            if ic_no_bed.sum() > 0:
                dead = self.ic_capacity * self.death_rate
                dead += ic_no_bed * self.ic_no_bed_death_rate
        self.dead += dead
        self.hospital_ic -= dead

        # hospital dead
        hospital_dead = hospitalized * self.hospital_death_rate
        self.dead += hospital_dead
        self.hospital -= hospital_dead

    def run(self, days):
        for day in range(days):
            # Do not use new data for calculating transitions
            infected = self.infected
            hospital = self.hospital
            exposed = self.exposed
            ic = self.hospital_ic

            # Transitions between compartments
            if self.scenario == 'party':
                if day < self.measure_day - 1:
                    self.infect(self.susceptible, infected, np.full((86, 86), 1))
                elif self.measure_day - 1 <= day < self.measure_day:
                    self.infect(self.susceptible, infected, np.array(
                        [[2 if 17 < i < 26 and 17 < j < 26 else 1 for i in range(86)] for j in range(86)]))
                elif self.measure_day <= day <= self.measure_day + 1:
                    self.infect(self.susceptible, infected, np.array(
                        [[2 * self.measure_factor if 17 < i < 26 and 17 < j < 26 else self.measure_factor for i in
                          range(86)] for j in range(86)]))
                else:
                    self.infect(self.susceptible, infected, np.full((86, 86), self.measure_factor))

            else:
                if day >= self.measure_day:
                    self.infect(self.susceptible, infected, np.full((86, 86), self.measure_factor))
                else:
                    self.infect(self.susceptible, infected, np.full((86, 86), 1))
            # self.infect(self.susceptible_hospital_staff, infected, 5)
            self.exp_to_inf(exposed)
            self.recover(infected, hospital, ic)
            self.go_to_hospital(infected)
            self.go_to_ic(hospital)
            self.die(hospital, ic)

            # Save data for graphing
            self.infected_data = np.append(self.infected_data, self.infected.sum())
            self.exposed_data = np.append(self.exposed_data, self.exposed.sum())
            self.recovered_data = np.append(self.recovered_data, self.recovered.sum())
            self.hospital_data = np.append(self.hospital_data, self.hospital.sum())
            self.ic_data = np.append(self.ic_data, self.hospital_ic.sum())
            self.dead_data = np.append(self.dead_data, self.dead.sum())
            self.case_data = np.append(self.case_data, self.infected.sum() + self.recovered.sum() + self.hospital.sum()
                                       + self.hospital_ic.sum() + self.dead.sum())
            self.hospital_total_data = np.append(self.hospital_total_data, self.hospital_total.sum())

            # print(self.susceptible.sum() + infected.sum() + self.recovered.sum())
            # print(i)
        # print('---')
        # print('Number of susceptible: ', self.susceptible.sum())
        # print('Number of infected: ', self.infected.sum())
        # print('Number of recovered: ', self.recovered.sum())

        # print(self.susceptible.sum() + self.infected.sum() + self.recovered.sum())
