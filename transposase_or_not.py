import sys

def read_transfile(filename):
    coords = []
    with open(filename, 'r') as f:
        for line in f:
            tmp = line.split()
            name = tmp[0]
            coord1 = int(tmp[1])
            coord2 = int(tmp[2])
            is_transposase = None
            if tmp[3] == "Non-transposase":
                is_transposase = False
            elif tmp[3] == "Transposase":
                is_transposase = True

            if (coord1 < coord2):
                coords.append((name, coord1, coord2, is_transposase))
            else:
                coords.append((name, coord1, coord2, is_transposase))
    return coords


def read_repeatfile(filename):
    coords = []
    with open(filename, 'r') as f:
        for line in f:
            tmp = line.split()
            coord1 = int(tmp[1])
            coord2 = int(tmp[2])
            length = int(tmp[3])
            if coord1 < coord2:
                coords.append((coord1, coord2, length))
            else:
                coords.append((coord1, coord2, length))
    return coords


def count_genes_between(repeat_coord1, repeat_coord2, genes_coords):
    count = 0
    count_of_transposase = 0
    genes_ids = []
    for (id, coord1, coord2, is_transposase) in genes_coords:
        if (coord1 > repeat_coord1) and (coord2 < repeat_coord2):
            count += 1
            if is_transposase:
                count_of_transposase += 1
            else:
                genes_ids.append(id)
    return count, count_of_transposase, genes_ids


def run(gene_filename, transposon_filename):
    genes_coords = read_transfile(gene_filename)
    repeats_coords = read_repeatfile(transposon_filename)

    information_file = open(gene_filename + "_Information.csv", 'w')
    information_file.write('repeat_start' + '\trepeat_end' + '\trepeat_length' + '\t#genes' + '\t#transposase' + '\tgenes_ids\n')


    with open("Count_between_repeats1.txt", 'w') as w:
        for i in range(len(repeats_coords)):
            (r1_coord1, r1_coord2, repeat_length) = repeats_coords[i]
            count, transp_count, genes_ids = count_genes_between(r1_coord1, r1_coord2, genes_coords)
            w.write('count of genes between ' + str(count) + ' count of transposase between ' + str(transp_count) + "\n")
            information_file.write(str(r1_coord1) + '\t' + str(r1_coord2) + '\t' + str(repeat_length) + '\t'
                                   + str(count) + '\t' + str(transp_count) + '\t')

            for gene_id in genes_ids:
                information_file.write(gene_id + ',')
            information_file.write('\n')

    information_file.close()

if __name__ == "__main__":
    gene_filename = sys.argv[1]
    repeat_filename = sys.argv[2]
    run(gene_filename, repeat_filename)
