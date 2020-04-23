from src.CCMatrix import CCMatrix
import numpy as np

# Eerste infected persoon beinvloed wss de curve

infectious_rate = 0.05
# incubation_rate = 1.0/7.0   # incubation period of 7 days
recovery_rate = 1.0/21.0    # infectious period of 21 days
contact_reduction = 5       # contact reduction after 30 days (x times less contacts)


def main():
    cc = CCMatrix('cc15.csv', 'leeftijden.csv')
    contact_matrix = cc.cc_matrix
    susceptible = cc.age_array
    exposed = np.zeros(86)
    infected = np.zeros(86)
    infected[35] = 2
    recovered = np.zeros(86)
    contact_matrix = contact_matrix * infectious_rate

    for i in range(60):
        x = np.transpose(np.asmatrix(susceptible))
        # susceptible_susceptible = np.matmul(x, np.asmatrix(susceptible))
        if i == 30:
            contact_matrix *= 1/contact_reduction
        infected_new, susceptible = infect(contact_matrix, infected.copy(), susceptible, x)
        recoveries = infected * recovery_rate
        recovered += recoveries
        infected_new -= recoveries
        infected_new = np.maximum(infected_new, np.zeros(86))

        infected = infected_new
        print(susceptible.sum() + infected.sum() + recovered.sum())
        print(i)

    print('---')
    print(susceptible.sum())
    print(infected.sum())
    print(recovered.sum())

    print(susceptible.sum() + infected.sum() + recovered.sum())


def infect(contact_matrix, infected, susceptible, x):
    susceptible_infected = np.matmul(x, np.asmatrix(infected))
    contacts = np.multiply(contact_matrix, susceptible_infected)
    contacts = np.asarray(contacts)

    for column in contacts.transpose():
        column = np.minimum(susceptible, column)
        infected += column
        susceptible -= column
    susceptible = np.maximum(susceptible, np.zeros(86))
    return infected, susceptible


if __name__ == "__main__":
    main()
