import numpy as np


def get_ages(attributes_filename):
    with open(attributes_filename) as file:
        file.readline()
        return np.asarray([int(line.split()[3]) for line in file.readlines()], dtype=int)


def get_genders_idxs(attributes_filename):
    with open(attributes_filename) as file:
        file.readline()
        genders = [line.split()[2] for line in file.readlines()]
        gender_idxs = {'M': [], 'F': []}

        for idx, gender in enumerate(genders):
            if gender == 'M':
                gender_idxs['M'].append(idx)
            elif gender == 'F':
                gender_idxs['F'].append(idx)
            else:
                raise ValueError('Icorrect gender')

        return gender_idxs['M'], gender_idxs['F']


def get_betas(betas_filename):
    with open(betas_filename) as file:
        file.readline()
        flag = False
        names = []
        betas = []

        while not flag:
            line = file.readline().split()
            if not len(line):
                flag = True
                break
            names.append(line.pop(0))
            betas.append(np.asarray([float(beta) for beta in line], dtype=float))

        return np.asarray(betas), names
