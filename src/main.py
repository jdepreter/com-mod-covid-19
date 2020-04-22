from src.CCMatrix import CCMatrix
import numpy as np

# Eerste infected persoon beinvloed wss de curve

transmission_rate = 0.5
recovery_rate = 0.1


def main():
    cc = CCMatrix('../data/cc15.csv', '../data/leeftijden.csv')
    contact_matrix = cc.cc_matrix
    susceptible = cc.age_array
    exposed = np.zeros(86)
    infected = np.zeros(86)
    infected[22] = 1
    infected[21] = 1
    recovered = np.zeros(86)
    contact_matrix = contact_matrix * transmission_rate

    for i in range(5):
        x = np.transpose(np.asmatrix(susceptible))
        # susceptible_susceptible = np.matmul(x, np.asmatrix(susceptible))
        infected_new, susceptible = infect(contact_matrix, infected.copy(), susceptible, x)
        recoveries = infected * recovery_rate
        recovered += recoveries
        infected_new -= recoveries
        infected_new = np.maximum(infected_new, np.zeros(86))

        infected = infected_new

    print(infected.sum())
    print(recovered.sum())


def infect(contact_matrix, infected, susceptible, x):
    susceptible_infected = np.matmul(x, np.asmatrix(infected))
    contacts = np.multiply(contact_matrix, susceptible_infected)
    contacts = np.asarray(contacts)
    # print(contacts)
    for column in contacts.transpose():
        infected += column
        susceptible -= column
    susceptible = np.maximum(susceptible, np.zeros(86))
    return infected, susceptible


if __name__ == "__main__":
    main()
