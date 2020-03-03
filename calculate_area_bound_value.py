import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt

from data import get_ages, get_genders_idxs, get_betas, get_bad_cpgs, get_gene_dict
from polygon import get_polygon, get_polygons_areas, get_polygon_xy


def calculate_area_bound_value(age_idx, gender_idx, data_path, base_name):
    attributes_filename   = data_path+f'attributes {base_name}.txt'
    # betas_filename        = data_path+f'test.txt'
    betas_filename        = data_path+f'{base_name}_average_beta.txt'
    bad_cpgs_filename     = data_path+'bad_cpgs.txt'
    annotations_filename  = data_path+'annotations.txt'

    ages = get_ages(attributes_filename, age_idx)
    male_idxs, female_idxs = get_genders_idxs(attributes_filename, gender_idx)
    bad_cpgs = get_bad_cpgs(bad_cpgs_filename, annotations_filename)
    betas, cpgs_names = get_betas(betas_filename, bad_cpgs)

    male_ages   = np.asarray([ages[i] for i in male_idxs])
    female_ages = np.asarray([ages[i] for i in female_idxs])

    areas = []

    for cpg_idx, _ in enumerate(cpgs_names):
        male_betas   = np.asarray([betas[cpg_idx, i] for i in male_idxs])
        female_betas = np.asarray([betas[cpg_idx, i] for i in female_idxs])

        male_regression = LinearRegression().fit(male_ages.reshape(-1, 1), male_betas)
        male_slope = male_regression.coef_
        male_intercept = male_regression.intercept_

        female_regression = LinearRegression().fit(female_ages.reshape(-1, 1), female_betas)
        female_slope = female_regression.coef_
        female_intercept = female_regression.intercept_

        male_polygon   = get_polygon(male_betas,   male_ages,   male_slope,   male_intercept)
        female_polygon = get_polygon(female_betas, female_ages, female_slope, female_intercept)

        intersection_area, union_area = get_polygons_areas(male_polygon, female_polygon)
        areas.append(intersection_area/union_area)

    q = .25
    quantile = np.quantile(areas, q)
    percentile = np.percentile(areas, q)

    print(f'quantile: {quantile}')
    print(f'percentile: {percentile}')

    return np.asarray(areas, dtype=float)


def load_areas():
    file = np.load('./results/areas.npz')
    areas40279 = file['areas40279']
    areas87571 = file['areas87571']

    print(f'{len(areas40279)}, {len(areas87571)}')
    print(f'{np.min(areas40279)}, {np.min(areas87571)}')
    print(f'{np.max(areas40279)}, {np.max(areas87571)}')

    plt.hist(areas40279, bins = 20)
    plt.savefig(f'./results/areas_GSE40279.png')
    plt.clf()

    plt.hist(areas87571, bins = 20)
    plt.savefig(f'./results/areas_GSE87571.png')
    plt.clf()


def main():
    data_path = '../../DNA_Methylation/methylation_data/'

    gse40279 = 'GSE40279'
    gse87571 = 'GSE87571'

    print(gse40279)
    areas40279 = calculate_area_bound_value(2, 3, data_path, gse40279)

    print(gse87571)
    areas87571 = calculate_area_bound_value(3, 2, data_path, gse87571)

    np.savez('./results/areas.npz', areas40279=areas40279, areas87571=areas87571)

    load_areas()


if __name__ == '__main__':
    main()
