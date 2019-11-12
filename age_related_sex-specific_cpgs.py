import argparse
from os.path import isfile

from metrics import fill_table, save_table, read_table
from tables import get_common_table


parser = argparse.ArgumentParser()
parser.add_argument('-dp', dest='data_path', type=str, nargs='*')
parser.add_argument('-sp', dest='save_path', type=str, nargs='*')


def main():
    data_path = '../../DNA_Methylation/methylation_data/'
    save_path = './results/'

    gse40279 = 'GSE40279'
    gse87571 = 'GSE87571'

    args = parser.parse_args()
    if args.data_path:
        data_path = args.data_path[0]
    if args.save_path:
        save_path = args.save_path[0]

    # GSE40279
    gse40279_table = None
    if isfile(save_path+gse40279+'_result_table.txt'):
        gse40279_table = read_table(save_path+gse40279+'_result_table.txt')
    else:
        gse40279_table = fill_table(2, 3, data_path, gse40279, save_path)

    print(f'Best CpGs ({gse40279}):')
    for j in range(min(10, len(gse40279_table))):
        print(gse40279_table[j])

    # GSE87571
    gse87571_table = None
    if isfile(save_path+gse87571+'_result_table.txt'):
        gse87571_table = read_table(save_path+gse87571+'_result_table.txt')
    else:
        gse87571_table = fill_table(3, 2, data_path, gse87571, save_path)

    print(f'Best CpGs ({gse87571}):')
    for j in range(min(10, len(gse87571_table))):
        print(gse87571_table[j])

    # Common CpGs
    common_table = get_common_table(gse40279_table, gse87571_table)

    print('Best common CpGs:')
    for i in range(min(10, len(common_table))):
        print(common_table[i])

if __name__ == '__main__':
    main()
