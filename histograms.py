import numpy as np
from matplotlib import pyplot as plt

from data import get_ages, get_genders_idxs


def get_histogram(age_idx, gender_idx, data_path, base_name, save_path, swap_flag=True):
    attributes_filename = data_path+f'{base_name}/attributes.txt'

    male_idxs, female_idxs = get_genders_idxs(attributes_filename, gender_idx)
    ages = get_ages(attributes_filename, age_idx)

    male_ages   = np.asarray([ages[i] for i in male_idxs])
    female_ages = np.asarray([ages[i] for i in female_idxs])

    if swap_flag:
        plt.hist(female_ages, bins=35, range=(0, 101), color='#FF69B4', label='gender(F)', alpha=.6)
        plt.hist(male_ages,   bins=35, range=(0, 101), color='#87CEFA', label='gender(M)', alpha=.6)
    else:
        plt.hist(male_ages,   bins=35, range=(0, 101), color='#87CEFA', label='gender(M)', alpha=.6)
        plt.hist(female_ages, bins=35, range=(0, 101), color='#FF69B4', label='gender(F)', alpha=.6)

    plt.title(f'{base_name}')
    plt.xlabel('age')
    plt.ylabel('count')
    plt.legend()

    plt.savefig(f'{save_path}histogram_{base_name}.png')

    plt.clf()
