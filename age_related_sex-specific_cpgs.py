# from matplotlib import pyplot as plt
# from shapely.geometry.polygon import Polygon


# def main():
#     left = Polygon(((0, 0), (2, 0), (2, 1), (0, 1)))
#     lx, ly = left.exterior.xy
#     right = Polygon(((1, 0), (3, 0), (3, 1), (1, 1)))
#     rx, ry = right.exterior.xy

#     intersection = left.intersection(right)
#     ix, iy = intersection.exterior.xy

#     print(f'areas: {left.area}, {right.area}, {intersection.area}')

#     plt.plot(lx, ly)
#     plt.plot(rx, ry)
#     plt.plot(ix, iy)
#     plt.show()


import argparse

from metrics import fill_table, save_table


parser = argparse.ArgumentParser()
parser.add_argument('-dp', dest='data_path', type=str, nargs='*')
parser.add_argument('-sp', dest='save_path', type=str, nargs='*')


def main():
    data_path = 'D:/Projects/DNA_Methylation/methylation_data/'
    save_path = 'D:/Projects/Python/Age_Related_Sex-Specific_CpGs/results/'

    args = parser.parse_args()
    if args.data_path:
        data_path = args.data_path[0]
    if args.save_path:
        save_path = args.save_path[0]

    metrics_table = fill_table(data_path+'attributes GSE87571.txt', data_path+'test.txt', save_path)
    metrics_table = fill_table(data_path+'attributes GSE87571.txt', data_path+'GSE87571_average_beta.txt', save_path)
    for i, row in enumerate(metrics_table):
        print(row)


if __name__ == '__main__':
    main()
