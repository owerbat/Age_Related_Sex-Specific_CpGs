def get_common_table(table1, table2):
    names1 = set([row[0] for row in table1])
    names2 = set([row[0] for row in table2])

    common_names = names1.intersection(names2)

    min_table = table1 if len(table1) < len(table2) else table2

    common_table = []
    for row in min_table:
        if row[0] in common_names:
            common_table.append(row)
    
    return common_table
