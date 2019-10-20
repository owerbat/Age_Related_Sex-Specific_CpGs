import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt

from data import get_ages, get_genders_idxs, get_betas
from polygon import get_polygon, get_polygons_areas, get_polygon_xy


def fill_table(attributes_filename, betas_filename, save_path):
    ages = get_ages(attributes_filename)
    male_idxs, female_idxs = get_genders_idxs(attributes_filename)
    betas, cpgs_names = get_betas(betas_filename)

    male_ages   = np.asarray([ages[i] for i in male_idxs])
    female_ages = np.asarray([ages[i] for i in female_idxs])

    table = []

    for cpg_idx, cpg_name in enumerate(cpgs_names):
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
        coef = intersection_area/union_area

        if (coef < .5) and ((np.abs(male_slope) > .001) or (np.abs(female_slope) > .001)):
            table.append((cpg_name, coef, male_slope[0], female_slope[0]))

            # observations
            plt.plot(male_ages,   male_betas, 'o', color='#87CEFA')
            plt.plot(female_ages, female_betas, 'o', color='#FF69B4')

            # linear regression lines
            plt.plot((np.min(male_ages), np.max(male_ages)),
                     (male_slope*np.min(male_ages)+male_intercept, male_slope*np.max(male_ages)+male_intercept),
                     color='#000080')
            plt.plot((np.min(female_ages), np.max(female_ages)),
                     (female_slope*np.min(female_ages)+female_intercept, female_slope*np.max(female_ages)+female_intercept),
                     color='#D71274')

            # polygons
            plt.plot(*get_polygon_xy(male_polygon), color='#0000FF')
            plt.plot(*get_polygon_xy(female_polygon), color='#FC0043')

            plt.title(cpg_name)
            plt.xlabel('age')
            plt.ylabel('beta')

            plt.savefig(save_path+'plots/'+cpg_name+'.png')

            plt.clf()

    table.sort(key=lambda items: items[1])

    save_table(table, save_path+'result_table.txt')

    return table


def save_table(table, save_filename):
    with open(save_filename, 'w') as file:
        for i, row in enumerate(table):
            file.write(f'{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n')
