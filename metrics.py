import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt

from data import get_ages, get_genders_idxs, get_betas, get_bad_cpgs, get_gene_dict
from polygon import get_polygon, get_polygons_areas, get_polygon_xy


# Intersection area bound values (25%)
AREA_BOUND_40279 = 0.8137639994808806
AREA_BOUND_87571 = 0.7914361897297456
AREA_BOUND_EPIC  = 0.7152952487992403
AREA_BOUND_55763 = 0.8820720148976503

area_bound_dict = {
    'GSE40279' : AREA_BOUND_40279,
    'GSE87571' : AREA_BOUND_87571,
    'epic'     : AREA_BOUND_EPIC,
    'GSE55763' : AREA_BOUND_55763
}


def fill_table(age_idx, gender_idx, data_path, base_name, save_path):
    attributes_filename   = data_path+f'{base_name}/attributes.txt'
    betas_filename        = data_path+f'{base_name}/betas.txt'
    bad_cpgs_filename     = data_path+f'{base_name}/bad_cpgs.txt'
    annotations_filename  = data_path+'annotations.txt'

    ages = get_ages(attributes_filename, age_idx)
    male_idxs, female_idxs = get_genders_idxs(attributes_filename, gender_idx)
    bad_cpgs = get_bad_cpgs(bad_cpgs_filename, annotations_filename)
    betas, cpgs_names = get_betas(betas_filename, bad_cpgs)
    gene_dict = get_gene_dict(annotations_filename)

    male_ages   = np.asarray([ages[i] for i in male_idxs])
    female_ages = np.asarray([ages[i] for i in female_idxs])

    table = []

    for cpg_idx, cpg_name in enumerate(cpgs_names):
        gene_name = None
        try:
            gene_name = gene_dict[cpg_name]
        except:
            gene_name = 'na'

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
        area = intersection_area/union_area

        area_bound = get_area_bound(base_name)

        if (area < area_bound) and ((np.abs(male_slope) > .001) or (np.abs(female_slope) > .001)):
            table.append((cpg_name, gene_name, area, male_slope[0], female_slope[0]))

            # observations
            plt.plot(male_ages,   male_betas,   'o', color='#87CEFA')
            plt.plot(female_ages, female_betas, 'o', color='#FF69B4')

            # linear regression lines
            plt.plot((np.min(male_ages), np.max(male_ages)),
                     (male_slope*np.min(male_ages)+male_intercept, male_slope*np.max(male_ages)+male_intercept),
                     color='#000080')
            plt.plot((np.min(female_ages), np.max(female_ages)),
                     (female_slope*np.min(female_ages)+female_intercept, female_slope*np.max(female_ages)+female_intercept),
                     color='#D71274')

            # polygons
            plt.plot(*get_polygon_xy(male_polygon),   color='#0000FF')
            plt.plot(*get_polygon_xy(female_polygon), color='#FC0043')

            plt.title(f'{cpg_name} ({gene_name})')
            plt.xlabel('age')
            plt.ylabel('beta')
            plt.axis((0, 110, -.1, 1.1))

            plt.savefig(f'{save_path}plots/{base_name}/{cpg_name}.png')

            plt.clf()

    table.sort(key=lambda items: items[2])

    save_table(table, save_path+f'{base_name}_result_table.txt')

    return table


def save_table(table, save_filename):
    with open(save_filename, 'w') as file:
        file.write('cpg\tgene\tintersection_area\tslope_male\tslope_female\n')
        for row in table:
            file.write(f'{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\n')


def read_table(table_filename):
    with open(table_filename) as file:
        table = []
        file.readline()
        for line in file:
            line = line.split()
            table.append((line[0], line[1], float(line[2]), float(line[3]), float(line[4])))
        return table


def get_area_bound(base_name):
    return area_bound_dict[base_name]
