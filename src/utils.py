import csv
def dump_matrix_to_csv(matrix, filename, delimiter=','):
    dense_matrix = matrix.toarray()
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter)
        for row in dense_matrix:
            writer.writerow(row)
    print(f"Matrix dumped to {filename}")