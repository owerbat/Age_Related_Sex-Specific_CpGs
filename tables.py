from numpy import mean, max


def get_common_table(table1, table2):
    names1 = set([row[0] for row in table1])
    names2 = set([row[0] for row in table2])

    common_names = names1.intersection(names2)

    dict1 = dict([(row[0], (row[1], row[2], row[3], row[4])) for row in table1])
    dict2 = dict([(row[0], (row[1], row[2], row[3], row[4])) for row in table2])

    common_table = []

    for name in common_names:
        common_table.append((name, dict1[name][0], (dict1[name][1] + dict2[name][1])/2,
                             max(dict1[name][2], dict2[name][2]), max(dict1[name][3], dict2[name][3])))

    return sorted(common_table, key=lambda item: item[2])


def get_common_table_multiple(tables):
    names = [set([row[0] for row in table]) for table in tables]
    common_names = names[0]

    for i in range(1, len(names)):
        common_names = common_names.intersection(names[i])

    data = [dict([(row[0], (row[1], row[2], row[3], row[4])) for row in table]) for table in tables]

    common_table = []

    for name in common_names:
        gene          = data[0][name][0]
        areas         = mean([item[name][1] for item in data])
        male_slopes   = max([abs(item[name][2]) for item in data])
        female_slopes = max([abs(item[name][3]) for item in data])

        common_table.append([name, gene, areas, male_slopes, female_slopes])

    common_table.sort(key=lambda items: items[2])

    return common_table


def get_area_table(tables):
    names = [set([row[0] for row in table]) for table in tables]
    common_names = names[0]

    for i in range(1, len(names)):
        common_names = common_names.intersection(names[i])

    data = [dict([(row[0], (row[1], row[2], row[3], row[4])) for row in table]) for table in tables]

    common_table = []

    for name in common_names:
        gene = data[0][name][0]
        res = [name, gene]

        for item in data:
            res.append(item[name][1])

        common_table.append(res)

    common_table.sort(key=lambda items: items[2])

    return common_table


def save_area_table(table, save_filename):
    with open(save_filename, 'w') as file:
        file.write('cpg\tgene\tGSE40279\tGSE87571\tEPIC\tGSE55763\n')
        for row in table:
            file.write(f'{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n')
