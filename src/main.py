from src.CCMatrix import CCMatrix
from src.Model import Model
from src.Grapher import Grapher
import matplotlib.pyplot as plt
import numpy as np

# Eerste infected persoon beinvloed wss de curve
# TODO: optimale infectious_rate bepalen a.d.h.v. optimal fit met data van Italie


def optimal_fit(reference_infected, infectious_rate, incubation_rate, recovery_rate, cc, num_days):
    # initial model, do once with higher infectious rate, once with lower, continue with best fit
    inf_rate_diff = 0.001
    model = Model(cc, infectious_rate + inf_rate_diff, incubation_rate, recovery_rate)
    model.run(num_days)
    residual = np.absolute(model.case_data - reference_infected)
    pos_sum_of_squares = np.sum(residual ** 2)

    inf_rate_diff = -0.001
    model = Model(cc, infectious_rate + inf_rate_diff, incubation_rate, recovery_rate)
    model.run(num_days)
    residual = np.absolute(model.case_data - reference_infected)
    neg_sum_of_squares = np.sum(residual ** 2)

    if neg_sum_of_squares < pos_sum_of_squares:
        infectious_rate -= 0.001
        prev_sum_of_squares = neg_sum_of_squares
        while True:
            infectious_rate -= 0.001
            print(infectious_rate)
            model = Model(cc, infectious_rate, incubation_rate, recovery_rate)
            model.run(num_days)
            residual = np.absolute(model.case_data - reference_infected)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares or infectious_rate <= 0 or infectious_rate >= 1:
                break
            else:
                prev_sum_of_squares = sum_of_squares
    else:
        infectious_rate += 0.001
        prev_sum_of_squares = pos_sum_of_squares
        while True:
            infectious_rate += 0.001
            print(infectious_rate)
            model = Model(cc, infectious_rate, incubation_rate, recovery_rate)
            model.run(num_days)
            residual = np.absolute(model.case_data - reference_infected)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares or infectious_rate <= 0 or infectious_rate >= 1:
                break
            else:
                prev_sum_of_squares = sum_of_squares
    return infectious_rate


def main():
    infectious_rate = 0.05
    incubation_rate = 1.0 / 3.0  # incubation period of 3 days
    recovery_rate = 1.0 / 6.0  # infectious period of 6 days
    cc = CCMatrix('cc15.csv', 'eurostat_pop_age.csv')

    reference_data = [0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,20,62,155,229,322,453,655,888,1128,
                       1694,2036,2502,3089,3858,4636,5883,7375,9172,10149,12462,12462,17660,21157,24747,27980,31506,
                       35713,41035,47021,53578,59138,63927,69176,74386,80589,86498,92472,97689,101739,105792,110574,
                       115242,119827,124632,128948,132547,135586,139422,143626,147577,152271,156363,159516,162488,165155,
                       168941,172434,175925,178972,181228,183957,187327,189973,192994,195351,197675,199414,201505,203591,
                       205463,207428,209328,210717][20:40]

    reference_cases = np.array(reference_data)

    infectious_rate = optimal_fit(reference_cases, infectious_rate, incubation_rate, recovery_rate, cc, 20)

    model = Model(cc, infectious_rate, incubation_rate, recovery_rate)
    days = 130
    model.run(days)

    # plt.plot(model.infected_data, label='Model infected')
    plt.plot(model.hospital_data, label='hospitalized')
    plt.plot(model.hospital_total_data, label='hospitalized total')
    plt.plot(model.ic_data, label='intensive care')
    # plt.plot(model.recovered_data, label='Model recovered')
    plt.plot(model.dead_data, label='died')
    # plt.plot(model.case_data, label='Model cases')
    # plt.plot(reference_cases, label='Reference cases')
    plt.legend()
    plt.show()
    print(model.dead_data[-1])
    # g = Grapher(days, [model.infected_data, model.hospital_data, model.ic_data, model.recovered_data, model.dead_data],
    #             ["Infected", "Hospital", "IC", "Recovered", "Dead"], display=True, save=True)
    # g = Grapher(days, [model.infected_data, model.hospital_data, model.ic_data, model.recovered_data, model.dead_data],
    #             ["Hospital", "IC", "Recovered", "Dead"], display=True, save=True)
    # g.animate("model")


if __name__ == "__main__":
    main()
