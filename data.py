import numpy as np


def get_ages(attributes_filename, age_idx):
    with open(attributes_filename) as file:
        file.readline()
        return np.asarray([int(line.split()[age_idx]) for line in file.readlines()], dtype=int)


def get_genders_idxs(attributes_filename, gender_idx):
    with open(attributes_filename) as file:
        file.readline()
        genders = [line.split()[gender_idx] for line in file.readlines()]
        gender_idxs = {'M': [], 'F': []}

        for idx, gender in enumerate(genders):
            if gender == 'M':
                gender_idxs['M'].append(idx)
            elif gender == 'F':
                gender_idxs['F'].append(idx)
            else:
                raise ValueError('Incorrect gender')

        return gender_idxs['M'], gender_idxs['F']


def get_betas(betas_filename, bad_cpgs):
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
            name = line.pop(0)
            if name not in bad_cpgs:
                names.append(name)
                betas.append(np.asarray([float(beta) for beta in line], dtype=float))

        return np.asarray(betas), names


def get_bad_cpgs(bad_cpgs_filename, annotations_filename):
    bad_cpgs = {}

    with open(bad_cpgs_filename) as bad_cpgs_file:
        for line in bad_cpgs_file:
            bad_cpgs.update({line.rstrip('\n'): 0})

    with open(annotations_filename) as annotations_file:
        for line in annotations_file:
            line = line.split('\t')
            name, map_info = line[0], line[2]
            if map_info in ('', 'X', 'Y', 'NA'):
                bad_cpgs.update({name: 0})

    return bad_cpgs


def get_gene_dict(annotations_filename):
    with open(annotations_filename) as file:
        gene_dict = {}
        for line in file:
            line = line.split('\t')
            cpg_name, gene_name = line[0], line[5].split(';')[0]
            if gene_name != '':
                gene_dict.update({cpg_name: gene_name})
        return gene_dict
