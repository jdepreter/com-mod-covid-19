from src.CCMatrix import CCMatrix
from src.Model import Model
import matplotlib.pyplot as plt
import numpy as np

# Eerste infected persoon beinvloed wss de curve


def main():
    infectious_rate = 0.05
    incubation_rate = 1.0 / 3.0  # incubation period of 7 days
    recovery_rate = 1.0 / 14.0  # infectious period of 21 days
    cc = CCMatrix('cc15.csv', 'leeftijden.csv')

    reference_infected = np.fromfunction(lambda x, y: 2**(x/2.34), (60, 1), dtype=int)

    # initial model, do once with higher infectious rate, once with lower, continue with best fit
    inf_rate_diff = 0.001
    model = Model(cc, infectious_rate + inf_rate_diff, incubation_rate, recovery_rate)
    model.run(60)
    residual = np.absolute(model.infected_data - reference_infected)
    pos_sum_of_squares = np.sum(residual**2)

    inf_rate_diff = -0.001
    model = Model(cc, infectious_rate + inf_rate_diff, incubation_rate, recovery_rate)
    model.run(60)
    residual = np.absolute(model.infected_data - reference_infected)
    neg_sum_of_squares = np.sum(residual ** 2)

    if neg_sum_of_squares < pos_sum_of_squares:
        infectious_rate -= 0.001
        prev_sum_of_squares = neg_sum_of_squares
        while True:
            infectious_rate -= 0.001
            print(infectious_rate)
            model = Model(cc, infectious_rate, incubation_rate, recovery_rate)
            model.run(60)
            residual = np.absolute(model.infected_data - reference_infected)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares or infectious_rate <= 0 or infectious_rate >= 1:
                break
    else:
        infectious_rate += 0.001
        prev_sum_of_squares = pos_sum_of_squares
        while True:
            infectious_rate += 0.001
            print(infectious_rate)
            model = Model(cc, infectious_rate, incubation_rate, recovery_rate)
            model.run(60)
            residual = np.absolute(model.infected_data - reference_infected)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares or infectious_rate <= 0 or infectious_rate >= 1:
                break


    plt.plot(model.infected_data, label='Model infected')
    plt.plot(model.recovered_data, label='Model recovered')
    plt.plot(reference_infected, label='Reference infected')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
