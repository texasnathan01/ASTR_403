import matplotlib.pyplot as plt
import numpy as np
import math as m

# General
c = 3e8 # Speed of light (m/s)

# Fill data structure with data from file

def get_data():
    data = []
    file = open('lcparam_full_long.txt', 'r')
    line = file.readline()
    while not (line == ''):
        line = file.readline()
        line = line[0:-1]
        parts = line.split(' ')
        data.append(parts)
    file.close()
    data = data[0:-1]
    return data

def filter_data(data):
    zcmb = []  # redshift
    mb = []  # observed B-band magnitude
    dmbb = []  # error of  observed B-band magnitude
    mu_LCDM = []  # theoretical distance modulus
    distances = []  # FIXME: May not be necessary
    velocities = []

    # Fill individual data columns from pulled data
    for x in data:
        zcmb.append(float(x[1]))
        mb.append(float(x[4]))
        dmbb.append(float(x[5]))  # FIXME: Column label not listed
        mu_LCDM.append(mb_to_mu_LCDM(float(x[4])))
        distances.append(mb_to_distance(float(x[4])))  # FIXME: May not be necessary
        velocities.append(zcmb_to_velocity(float(x[1])))

    return zcmb, mb, dmbb, mu_LCDM, distances, velocities

# 1.
#Calculates distance modulus from mb using distance modulus formula
def mb_to_mu_LCDM(mb):
    return mb + 19.3

#Calculates distance in MPC from mb using distance modulus formula
def mb_to_distance(mb):
    exponent = ((mb + 19.3) / 5) - 5
    return 10 ** exponent

#Calculates velocity in km/s
def zcmb_to_velocity(zcmb):
    # Non relativistic
    return zcmb * c / 1000

    # Relativistic
    z1 = (zcmb + 1) ** 2
    return (z1 - 1)/(z1 + 1) * c / 1000


# 2.
def get_q0(omega_m, omega_lambda):
    return (1/2) * (omega_m) - (omega_lambda)

def get_dl(H0, z, q0):
    return (c * z / H0) * (1 + z*(1 - q0)/2)

def get_distance_modulus(dl):
    return 5 * m.log10(dl) + 25

def get_distance_moduluses(H0, omega_m, omega_l, zcmb):
    distance_moduluses = []
    for z in zcmb:
        q0 = get_q0(omega_m, omega_l)
        dl = get_dl(H0, z, q0)
        distance_modulus = get_distance_modulus(dl)
        distance_moduluses.append(distance_modulus)
    return distance_moduluses

# 3.
def get_differences(mb, distance_moduluses):
    differences = []
    for i in range(len(mb)):
        mb_x = mb[i]
        distance_modulus_x = distance_moduluses[i]
        differences.append(mb_x - distance_modulus_x)
    return differences

# 4.
def chisquare(observed, expected, error):
    sum = 0
    for i in range(len(observed)):
        sum += ((observed[i] - expected[i]) ** 2) / (error[i])
    return sum

