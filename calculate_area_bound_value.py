import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from data import get_ages, get_genders_idxs, get_betas, get_bad_cpgs, get_gene_dict
from polygon import get_polygon, get_polygons_areas, get_polygon_xy


# Intersection area bound values (25%)
AREA_BOUND_40279 = 0.8137639994808806
AREA_BOUND_87571 = 0.7914361897297456
AREA_BOUND_EPIC  = 0.7152952487992403
AREA_BOUND_55763 = 0.8820720148976503


def calculate_area_bound_value(age_idx, gender_idx, data_path, base_name):
    attributes_filename   = data_path+f'{base_name}/attributes.txt'
    betas_filename        = data_path+f'{base_name}/betas.txt'
    bad_cpgs_filename     = data_path+f'{base_name}/bad_cpgs.txt'
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

        intersection_area, union_area, _ = get_polygons_areas(male_polygon, female_polygon)
        areas.append(intersection_area/union_area)

    q = .25
    quantile = np.quantile(areas, q)
    percentile = np.percentile(areas, q)

    print(f'quantile: {quantile}')
    print(f'percentile: {percentile}')

    return np.asarray(areas, dtype=float)


def get_plot(values, bound_value, base_name, file_name):
    a, b, _ = plt.hist(values, bins = 100, density=True, histtype='step', range=(0, 1))
    plt.clf()

    centers = [(b[i]+b[i+1])/2 for i in range(len(a))]
    f = interp1d(centers, a, kind='cubic', bounds_error=False, fill_value=(a[0], a[-1]))
    x_int = np.arange(0, 1, .001)
    y_int = f(x_int)

    plt.plot(x_int, y_int, color='green', label='PDF')
    plt.vlines(bound_value, ymin=np.min(y_int), ymax=np.max(y_int), colors='red', label='квартиль 25%')
    a, b, _ = plt.hist(values, bins = 20, density=True, histtype='step', range=(0, 1), label='гистограмма')

    plt.xlabel('площадь')
    plt.title(base_name)
    plt.legend()
    plt.savefig(file_name)
    plt.clf()


def load_areas():
    file = np.load('./results/areas.npz')
    areas40279 = file['areas40279']
    areas87571 = file['areas87571']
    areasEPIC  = file['areasEPIC']
    areas55763 = file['areas55763']

    print(f'{len(areas40279)}, {len(areas87571)}, {len(areasEPIC)}, {len(areas55763)}')
    print(f'{np.min(areas40279)}, {np.min(areas87571)}, {np.min(areasEPIC)}, {np.min(areas55763)}')
    print(f'{np.max(areas40279)}, {np.max(areas87571)}, {np.max(areasEPIC)}, {np.max(areas55763)}')

    get_plot(areas40279, AREA_BOUND_40279, 'GSE40279', './results/areas_GSE40279.png')
    get_plot(areas87571, AREA_BOUND_87571, 'GSE87571', './results/areas_GSE87571.png')
    get_plot(areasEPIC,  AREA_BOUND_EPIC,  'epic',     './results/areas_epic.png')
    get_plot(areas55763, AREA_BOUND_55763, 'GSE55763', './results/areas_GSE55763.png')


def main():
    data_path = '../../DNA_Methylation/methylation_data/'

    gse40279 = 'GSE40279'
    gse87571 = 'GSE87571'
    epic     = 'epic'
    gse55763 = 'GSE55763'

    print(gse40279)
    areas40279 = calculate_area_bound_value(2, 3, data_path, gse40279)

    print(gse87571)
    areas87571 = calculate_area_bound_value(3, 2, data_path, gse87571)

    print(epic)
    areasEPIC  = calculate_area_bound_value(2, 3, data_path, epic)

    print(gse55763)
    areas55763 = calculate_area_bound_value(3, 2, data_path, gse55763)

    np.savez('./results/areas.npz', areas40279=areas40279, areas87571=areas87571,
             areasEPIC=areasEPIC, areas55763=areas55763)

    load_areas()


if __name__ == '__main__':
    main()
