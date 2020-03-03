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
    names = [set([row[0] for row in table for table in tables])]
    common_names = names[0]

    for i in range(1, len(names)):
        common_names = common_names.intersection(names[i])

    data = [dict([(row[0], (row[1], row[2], row[3], row[4])) for row in table]) for table in tables]

    common_table = []

    for name in common_names:
        gene          = data[0][name][0]
        areas         = mean([item[name][1] for item in data])
        male_slopes   = max([item[name][2] for item in data])
        female_slopes = max([item[name][3] for item in data])

        common_table.append([name, gene, areas, male_slopes, female_slopes])

    common_table.sort(key=lambda items: items[2])

    return common_table


def get_common_table_multiple(tables):
    names = [set([row[0] for row in table for table in tables])]
    common_names = names[0]

    for i in range(1, len(names)):
        common_names = common_names.intersection(names[i])

    data = [{} for table in tables]

    for idx, table in enumerate(tables):
        for row in table:
            data[idx].update({row[0]: row[1:]})

    common_table = []

    for name in common_names:
        gene          = data[0][name][0]
        areas         = mean([item[name][1] for item in data])
        male_slopes   = mean([item[name][2] for item in data])
        female_slopes = mean([item[name][3] for item in data])

        common_table.append([name, gene, areas, male_slopes, female_slopes])

    common_table.sort(key=lambda items: items[2])

    return common_table
