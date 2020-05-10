from src.CCMatrix import CCMatrix
from src.Model import Model
from src.Grapher import plot, animate
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import dirname, abspath

img_folder = os.path.join(dirname(dirname(abspath(__file__))), 'img')

# Eerste infected persoon beinvloed wss de curve

reference_cases_italy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 20,
                         62,
                         155, 229, 322, 453, 655, 888, 1128,
                         1694, 2036, 2502, 3089, 3858, 4636, 5883, 7375, 9172, 10149, 12462, 12462, 17660, 21157, 24747,
                         27980, 31506,
                         35713, 41035, 47021, 53578, 59138, 63927, 69176, 74386, 80589, 86498, 92472, 97689, 101739,
                         105792, 110574,
                         115242, 119827, 124632, 128948, 132547, 135586, 139422, 143626, 147577, 152271, 156363, 159516,
                         162488, 165155,
                         168941, 172434, 175925, 178972, 181228, 183957, 187327, 189973, 192994, 195351, 197675, 199414,
                         201505, 203591,
                         205463, 207428, 209328, 210717]


def optimal_fit(reference_infected, infectious_rate, contact_matrix, susceptible, num_days):
    # initial model, do once with higher infectious rate, once with lower, continue with best fit
    inf_rate_diff = 0.001
    model = Model(contact_matrix, susceptible, infectious_rate + inf_rate_diff)
    model.run(num_days)
    residual = np.absolute(model.case_data - reference_infected)
    pos_sum_of_squares = np.sum(residual ** 2)

    inf_rate_diff = -0.001
    model = Model(contact_matrix, susceptible, infectious_rate + inf_rate_diff)
    model.run(num_days)
    residual = np.absolute(model.case_data - reference_infected)
    neg_sum_of_squares = np.sum(residual ** 2)

    if neg_sum_of_squares < pos_sum_of_squares:
        infectious_rate -= 0.001
        prev_sum_of_squares = neg_sum_of_squares
        while True:
            infectious_rate -= 0.001
            model = Model(contact_matrix, susceptible, infectious_rate)
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
            model = Model(contact_matrix, susceptible, infectious_rate)
            model.run(num_days)
            residual = np.absolute(model.case_data - reference_infected)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares or infectious_rate <= 0 or infectious_rate >= 1:
                break
            else:
                prev_sum_of_squares = sum_of_squares
    return infectious_rate


def find_infectious_rate(contact_matrix, susceptible, start_infectious_rate, reference_cases):
    infectious_rate = optimal_fit(reference_cases, start_infectious_rate, contact_matrix, susceptible, 20)

    # # Plot with this infectious rate
    # plot_model(cc, infectious_rate, 20, reference_cases)

    return infectious_rate


def find_offset(contact_matrix, susceptible, reference_hospital, infectious_rate, days, start_offset=20,
                measure_day=0, first_patient_age=38):
    prev = float('inf')
    offset = start_offset
    measure_day, factor, sum_of_squares = find_lockdown_factor(contact_matrix, susceptible, reference_hospital,
                                                               infectious_rate, days, offset - 1,
                                                               measure_day, first_patient_age)

    measure_day, factor, sum_of_squares2 = find_lockdown_factor(contact_matrix, susceptible, reference_hospital,
                                                                infectious_rate, days, offset + 1,
                                                                measure_day, first_patient_age)
    diff = 1
    if sum_of_squares < sum_of_squares2:
        diff = -1

    counter = 1
    while sum_of_squares <= prev:
        print(counter, sum_of_squares)
        prev = sum_of_squares
        measure_day, factor, sum_of_squares = find_lockdown_factor(contact_matrix, susceptible, reference_hospital,
                                                                   infectious_rate,
                                                                   days, offset + counter * diff,
                                                                   measure_day, first_patient_age)
        counter += 1
    print(factor, measure_day + (counter - 1) * diff, sum_of_squares)

    return offset + (counter - 1) * diff, measure_day, factor


def find_lockdown_factor(contact_matrix, susceptible, reference_hospital, infectious_rate, days, offset=20,
                         measure_day=0, first_patient_age=38):
    prev = float('inf')
    factor, sum_of_squares = calculate_sum_of_squares(contact_matrix, days, infectious_rate, measure_day - 1, offset,
                                                      reference_hospital,
                                                      susceptible, first_patient_age)

    factor, sum_of_squares2 = calculate_sum_of_squares(contact_matrix, days, infectious_rate, measure_day + 1, offset,
                                                       reference_hospital,
                                                       susceptible, first_patient_age)
    diff = 1
    if sum_of_squares < sum_of_squares2:
        diff = -1

    counter = 1
    while sum_of_squares <= prev:
        print(counter, sum_of_squares)
        prev = sum_of_squares
        factor, sum_of_squares = calculate_sum_of_squares(contact_matrix, days, infectious_rate,
                                                          measure_day + counter * diff, offset, reference_hospital,
                                                          susceptible, first_patient_age)
        counter += 1
    print("offset:", offset, factor, measure_day + (counter - 1) * diff, sum_of_squares)

    return measure_day + (counter - 1) * diff, factor, sum_of_squares


def calculate_sum_of_squares(contact_matrix, days, infectious_rate, measure_day, offset, reference_hospital,
                             susceptible, first_patient_age):
    num_days = days + offset
    factor = 0.2
    factor_diff = 0.001
    model = Model(contact_matrix, susceptible, infectious_rate, measure_factor=factor + factor_diff,
                  measure_day=measure_day, first_patient_age=first_patient_age)
    model.run(num_days)
    temp = model.hospital_data + model.ic_data
    residual = np.absolute(temp[offset:] - reference_hospital)
    pos_sum_of_squares = np.sum(residual ** 2)
    factor_diff = -0.001
    model = Model(contact_matrix, susceptible, infectious_rate, measure_factor=factor + factor_diff,
                  measure_day=measure_day, first_patient_age=first_patient_age)
    model.run(num_days)
    temp = model.hospital_data + model.ic_data
    residual = np.absolute(temp[offset:] - reference_hospital)
    neg_sum_of_squares = np.sum(residual ** 2)
    if neg_sum_of_squares < pos_sum_of_squares:
        factor_diff -= 0.001
        prev_sum_of_squares = neg_sum_of_squares
        while True:
            factor_diff -= 0.001
            model = Model(contact_matrix, susceptible, infectious_rate, measure_factor=factor + factor_diff,
                          measure_day=measure_day, first_patient_age=first_patient_age)
            model.run(num_days)
            temp = model.hospital_data + model.ic_data
            residual = np.absolute(temp[offset:] - reference_hospital)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares or infectious_rate <= 0 or infectious_rate >= 1:
                break
            else:
                prev_sum_of_squares = sum_of_squares
    else:
        factor_diff += 0.001
        prev_sum_of_squares = pos_sum_of_squares
        while True:
            factor_diff += 0.001
            model = Model(contact_matrix, susceptible, infectious_rate, measure_factor=factor + factor_diff,
                          measure_day=measure_day, first_patient_age=first_patient_age)
            model.run(num_days)
            temp = model.hospital_data + model.ic_data
            residual = np.absolute(temp[offset:] - reference_hospital)
            sum_of_squares = np.sum(residual ** 2)
            if sum_of_squares > prev_sum_of_squares:
                break
            else:
                prev_sum_of_squares = sum_of_squares
    return factor + factor_diff, prev_sum_of_squares


def plot_model(contact_matrix, susceptible, infectious_rate, days, reference_cases=None, measure_factor: float = 1,
               measure_day=0, reference_hospital=None, offset=0, first_patient_age=38, name='model',
               scenario=None, cap_ic=False) -> Model:
    model = Model(contact_matrix, susceptible, infectious_rate, measure_factor=measure_factor, measure_day=measure_day,
                  first_patient_age=first_patient_age, scenario=scenario, cap_ic=cap_ic)
    model.run(days)

    y = [model.infected_data, model.recovered_data, model.exposed_data, model.hospital_data, model.ic_data,
         model.dead_data]
    labels = ['infected', 'recovered', 'exposed', 'hospitalized', 'ic', 'dead']
    plot(y, labels, name, 'Amount')

    if reference_cases is not None:
        y = [reference_cases, model.case_data]
        labels = ['reference cases', 'model cases']
        plot(y, labels, 'ref_cases', 'Amount')

    if reference_hospital is not None:
        y = [np.append(np.zeros(offset), reference_hospital), model.hospital_data + model.ic_data]
        labels = ['reference hospital + ic', 'model hospital + ic']
        plot(y, labels, 'ref_hospital', 'Amount')

    return model


def main():
    infectious_rate = .05
    days = 150

    cc = CCMatrix('cc15.csv', 'eurostat_pop_age.csv', 'COVID19BE_HOSP.csv')
    contact_matrix = cc.cc_matrix

    # Data van Italië, we gebruiken hier dag 20 tot 40 (van dag waarop aantal gevallen begint te stijgen tot dag waarop maatregelen genomen zijn)
    # reference_cases = np.array(reference_cases_italy[20:40])
    # reference_cases = reference_cases/(cc.italy_count.sum()/cc.belgium_count.sum())  # rescale to belgium population

    # Find 'best fit' infectious rate and plot a run of the model with this rate
    # infectious_rate = find_infectious_rate(cc.cc_matrix, cc.italic_belgium_count, infectious_rate, reference_cases)
    infectious_rate = 0.109
    print('Best fit infectious rate:', infectious_rate)

    # Plot with this rate
    # model = plot_model(contact_matrix, cc.belgium_count, infectious_rate=infectious_rate, days=21,
    #                    reference_cases=reference_cases, first_patient_age=45, name='test')

    # measure_day = 37
    # offset, measure_day, factor = find_offset(contact_matrix, cc.belgium_count, cc.belgium_hospital, infectious_rate,
    #                               len(cc.belgium_hospital), measure_day=measure_day, start_offset=20,
    #                                           first_patient_age=45)
    #
    # print("uitkomst", offset, measure_day, factor)
    offset = 25
    measure_day = 38
    factor = 0.1

    print('Best fit offset', offset, 'Day Measures Introduced', measure_day, 'Factor', factor)
    model = plot_model(contact_matrix, cc.belgium_count, infectious_rate=infectious_rate, days=days,
                       measure_factor=factor, measure_day=measure_day, reference_hospital=cc.belgium_hospital,
                       offset=offset, first_patient_age=45)
    print('Dead people:', model.dead_data[-1])
    print('Hospital people:', model.hospital_total_data[-1])
    print('IC people:', model.ic_data[-1])

    model_party = plot_model(contact_matrix, cc.belgium_count, infectious_rate=infectious_rate, days=days,
                             measure_factor=factor, measure_day=measure_day, reference_hospital=cc.belgium_hospital,
                             offset=offset, scenario='party', first_patient_age=45, name='model_party')
    print('Dead people:', model_party.dead_data[-1])
    print('Hospital people:', model_party.hospital_total_data[-1])
    print('IC people:', model_party.ic_data[-1])

    plot([model.hospital_data, model_party.hospital_data], ['Normal', 'Lockdown Parties'], 'hospital_diff',
         'Hospitalizations')
    plot([model.ic_data, model_party.ic_data], ['Normal', 'Lockdown Parties'], 'ic_diff', 'Intensive Care')
    plot([model.infected_data, model_party.infected_data], ['Normal', 'Lockdown Parties'], 'infected_diff',
         'Infections')
    plot([model.dead_data, model_party.dead_data], ['Normal', 'Lockdown Parties'], 'death_diff', 'Dead')

    plot([np.append(np.zeros(offset), cc.belgium_hospital), model_party.hospital_data + model_party.ic_data],
         ["Hospital + IC (parties)", "Reference Hospital + IC"], 'party_hospital_ic_diff', 'Hospitalized')

    # Uncomment for plots
    model_no_lockdown = plot_model(contact_matrix, cc.belgium_count, infectious_rate=infectious_rate, days=days,
                                   first_patient_age=45, cap_ic=False)

    plot([model_no_lockdown.infected_data, model_no_lockdown.hospital_data, model_no_lockdown.ic_data,
          model_no_lockdown.dead_data], ["Infected", "Hospital", "IC", "Dead"],
         name='model_no_lockdown', y_label='Amount')

    plot([model_no_lockdown.infected_data, model.infected_data], ["No lockdown", "Lockdown"],
         name='Infected_lockdown_diff', y_label='Infections')

    plot([model_no_lockdown.hospital_data, model.hospital_data], ["No lockdown", "Lockdown"],
         name='Hospital_lockdown_diff', y_label='Hospitalizations')

    plot([model_no_lockdown.dead_data, model.dead_data], ["No lockdown", "Lockdown"],
         name='Dead_lockdown_diff', y_label='Dead')

    plot([model_no_lockdown.ic_data, model.ic_data], ["No lockdown", "Lockdown"],
         name='IC_lockdown_diff', y_label='Intensive Care')

    print('Dead people:', model_no_lockdown.dead_data[-1])
    print('Hospital people:', model_no_lockdown.hospital_total_data[-1])
    print('IC people:', model_no_lockdown.ic_data[-1])

    # Uncomment for plots (icu cap)
    model_no_lockdown_cap_ic = plot_model(contact_matrix, cc.belgium_count, infectious_rate=infectious_rate, days=days,
                                   first_patient_age=45, cap_ic=True)

    plot([model_no_lockdown_cap_ic.infected_data, model_no_lockdown_cap_ic.hospital_data, model_no_lockdown_cap_ic.ic_data,
          model_no_lockdown_cap_ic.dead_data, model_no_lockdown_cap_ic.ic_overflow_data], ["Infected", "Hospital", "IC", "Dead", "IC Overflow"],
         name='model_no_lockdown-cap-ic', y_label='Amount')

    # plot([model_no_lockdown.ic_data, model_no_lockdown.ic_overflow_data],
    #      ["IC", "IC Overflow"],
    #      name='model_no_lockdown-cap-ic', y_label='Amount')

    plot([model_no_lockdown_cap_ic.infected_data, model_no_lockdown.infected_data], ["No lockdown (IC overflow)", "No lockdown"],
         name='Infected_lockdown_diff-cap-ic', y_label='Infections')

    plot([model_no_lockdown_cap_ic.hospital_data, model_no_lockdown.hospital_data], ["No lockdown (IC overflow)", "No lockdown"],
         name='Hospital_lockdown_diff-cap-ic', y_label='Hospitalizations')

    plot([model_no_lockdown_cap_ic.dead_data, model_no_lockdown.dead_data], ["No lockdown (IC overflow)", "No lockdown"],
         name='Dead_lockdown_diff-cap-ic', y_label='Dead')

    plot([model_no_lockdown_cap_ic.ic_data, model_no_lockdown.ic_data], ["No lockdown (IC overflow)", "No lockdown"],
         name='IC_lockdown_diff-cap-ic', y_label='Intensive Care')

    print('Dead people:', model_no_lockdown_cap_ic.dead_data[-1])
    print('Hospital people:', model_no_lockdown_cap_ic.hospital_total_data[-1])
    print('IC people:', model_no_lockdown_cap_ic.ic_data[-1])

    # gifs
    # animate([model.infected_data, model.hospital_data, model.ic_data, model.dead_data],
    #         ["Infected", "Hospital", "IC", "Dead"], display=False, save=True, name='model', y_label='Amount')
    #
    # animate([model.infected_data, model.hospital_data, model.ic_data, model.dead_data],
    #         ["Infected", "Hospital", "IC", "Dead"], display=False, save=True, name='model', y_label='Amount')
    #
    # animate([model_no_lockdown.infected_data, model_no_lockdown.hospital_data, model_no_lockdown.ic_data,
    #          model_no_lockdown.dead_data], ["Infected", "Hospital", "IC", "Dead"], display=False, save=True,
    #         name='model_no_lockdown', y_label='Amount')
    #
    # animate([model_no_lockdown.infected_data, model.infected_data], ["No lockdown", "Lockdown"],
    #         display=False, save=True, name='Infected_lockdown_diff', y_label='Infections')
    #
    # animate([model_no_lockdown.hospital_data, model.hospital_data], ["No lockdown", "Lockdown"],
    #         display=False, save=True, name='Hospital_lockdown_diff', y_label='Hospitalizations')
    #
    # animate([model_no_lockdown.dead_data, model.dead_data], ["No lockdown", "Lockdown"],
    #         display=False, save=True, name='Dead_lockdown_diff', y_label='Dead')
    #
    # animate([model_no_lockdown.ic_data, model.ic_data], ["No lockdown", "Lockdown"],
    #         display=False, save=True, name='IC_lockdown_diff', y_label='Intensive Care')


if __name__ == "__main__":
    main()
