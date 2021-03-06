import argparse
from os.path import isfile

from metrics import fill_table, save_table, read_table
from tables import get_common_table_multiple, get_area_table, save_area_table
from histograms import get_histogram


parser = argparse.ArgumentParser()
parser.add_argument('-dp', dest='DATA_PATH', type=str, nargs='*')
parser.add_argument('-sp', dest='SAVE_PATH', type=str, nargs='*')


DATA_PATH = '../../DNA_Methylation/methylation_data/'
SAVE_PATH = './results/'


def compute(base_name, age_idx, gender_idx):
    table = None
    if isfile(SAVE_PATH+base_name+'_result_table.txt'):
        table = read_table(SAVE_PATH+base_name+'_result_table.txt')
    else:
        table = fill_table(age_idx, gender_idx, DATA_PATH, base_name, SAVE_PATH)

    print(f'Best CpGs ({base_name}):')
    for j in range(min(10, len(table))):
        print(table[j])

    return table


def main():
    global DATA_PATH
    global SAVE_PATH

    args = parser.parse_args()
    if args.DATA_PATH:
        DATA_PATH = args.DATA_PATH[0]
    if args.SAVE_PATH:
        SAVE_PATH = args.SAVE_PATH[0]

    gse40279 = 'GSE40279'
    gse87571 = 'GSE87571'
    epic     = 'epic'
    gse55763 = 'GSE55763'
    results  = []

    # get_histogram(2, 3, DATA_PATH, gse40279, SAVE_PATH)
    # get_histogram(3, 2, DATA_PATH, gse87571, SAVE_PATH)
    # get_histogram(2, 3, DATA_PATH, epic,     SAVE_PATH)
    # get_histogram(3, 2, DATA_PATH, gse55763, SAVE_PATH, False)

    results.append(compute(gse40279, 2, 3))
    results.append(compute(gse87571, 3, 2))
    results.append(compute(epic,     2, 3))
    results.append(compute(gse55763, 3, 2))

    # Common CpGs
    common_table = get_common_table_multiple(results)
    save_table(common_table, SAVE_PATH+'common_result_table.txt')

    area_table = get_area_table(results)
    save_area_table(area_table, SAVE_PATH+'area_result_table.txt')

    print('Best common CpGs:')
    for i in range(min(10, len(common_table))):
        print(common_table[i])

if __name__ == '__main__':
    main()
