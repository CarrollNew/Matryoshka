import time
from itertools import zip_longest, islice
import sys
import matplotlib.pyplot as plt
import numpy as np


def genome_read(filename):
    genome_seq  = ''
    genome_name = ''
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line[0] == '>':
                genome_name = line[:-1]
            else:
                genome_seq += line[:-1]

    return genome_name, genome_seq


# можно ускорить
def to_int_keys(l):
    """
    l: iterable of keys
    returns: a list with integer keys
    """
    index = {v: i for i, v in enumerate(sorted(set(l)))}
    return [index[v] for v in l]


def suffix_matrix_smart(to_int_keys, s):
    n = len(s)
    k = 1
    line = to_int_keys(s)
    ans = [line]
    while max(line) < n - 1:
        line = to_int_keys(
            list(zip_longest(line, islice(line, k, None),
                             fillvalue=-1)))
        ans.append(line)
        k <<= 1
    return ans


def lcp(sm, i, j):
    """
    longest common prefix
    O(log(n))

    sm: suffix matrix
    """
    n = len(sm[-1])
    if i == j:
        return n - i
    k = 1 << (len(sm) - 2)
    ans = 0
    for line in sm[-2::-1]:
        if i >= n or j >= n:
            break
        if line[i] == line[j]:
            ans ^= k
            i += k
            j += k
        k >>= 1
    return ans

def create_all_3mers(array_of_3mers, s=''):
    if len(s) == 3:
        array_of_3mers.append(s)
        return
    alphabet = "ACGT"
    for letter in alphabet:
        create_all_3mers(array_of_3mers, s + letter)



# в том числе вернем пары массивов троек (пара ([1,7,8,9], [2,3,4,5,6,10])
# то есть сначала идут индексы левой стороны, затем правой
def create_index_massive_list(suffix_array, joined_string, lower_bound):
    index_list = [0 for _ in range(len(suffix_array))]
    is_left    = [False for _ in range(len(suffix_array))]

    # запишем пару строка и индекс еев куске генома целом
    for i, index in enumerate(suffix_array):
        if i < lower_bound:
            index_list[index] = (joined_string[i:], i)
            is_left[index] = True
        else:
            index_list[index] = (joined_string[i:], i)


    # теперь мы знаем какие суффиксы в каком порядке
    # и знаем какой суффикс какой части принадлежит (левой или правой)
    # найдем разделители массива по тройкам букв

    pair = ([], [])
    pairs_array = []

    all_kmers = []
    create_all_3mers(all_kmers)
    d = {}
    for kmer in all_kmers:
        d[kmer] = ([], [])

    for i, (suffix, index_in_genome) in enumerate(index_list):
        if suffix[:3].find('#') == -1 and len(suffix) > 2:
                if is_left[i]:
                    d[suffix[:3]][0].append(index_in_genome)
                else:
                    d[suffix[:3]][1].append(index_in_genome)

    pairs_array.append(pair)

    return index_list, pairs_array, d


def create_suffix_matrix(genome, beg_pos, end_pos=None):
    if end_pos is None:
        end_pos = beg_pos

    genome_first_part  = genome[:beg_pos]
    genome_second_part = genome[end_pos:]

    # сколько берем слева от позиции
    lower_bound_left_part = 5000

    # как далеко заходим вправо от позиции
    upper_bound_right_part = 5000

    first_part_last150 = genome_first_part[-lower_bound_left_part:]
    second_part_first150 = genome_second_part[:upper_bound_right_part]

    joined_string = first_part_last150 + '#' + second_part_first150
    suff_matrix = suffix_matrix_smart(to_int_keys, joined_string)

    return joined_string, suff_matrix, lower_bound_left_part


def compliment(seq):
    compl_seq = ''
    for c in seq:
        if c == 'A':
            compl_seq += 'T'
        elif c == 'C':
            compl_seq += 'G'
        elif c == 'G':
            compl_seq += 'C'
        elif c == 'T':
            compl_seq += 'A'
        else:
            raise Exception()
    return compl_seq


def create_reverse_compl_suffix_matrix(genome, beg_pos, end_pos=None):
    if end_pos is None:
        end_pos = beg_pos

    genome_first_part  = genome[:beg_pos]
    genome_second_part = genome[end_pos:]

    # сколько берем слева от позиции
    lower_bound_left_part = 5000

    # как далеко заходим вправо от позиции
    upper_bound_right_part = 5000

    first_part_last150 = genome_first_part[-lower_bound_left_part:]
    second_part_first150 = compliment(reversed(genome_second_part[:upper_bound_right_part]))

    joined_string = first_part_last150 + '#' + second_part_first150
    suff_matrix = suffix_matrix_smart(to_int_keys, joined_string)

    return joined_string, suff_matrix, lower_bound_left_part, upper_bound_right_part

def find_direct_repeats(genome, beg_pos, end_pos):
    joined_string, suff_matrix, lower_bound_left_part = create_suffix_matrix(genome, beg_pos, end_pos)
    suff_array = suff_matrix[-1]

    print('BE', beg_pos, end_pos)

    # index_list - суффиксы в лексикограяичесокм порядке
    # pairs_array = массив пар, где в паре (массив с индексами слева от pos c 3меров в начале)
    # справа массив с индексами спрва от pos с тем же 3мером вначале
    index_list, pairs_array, d = create_index_massive_list(suff_array, joined_string, lower_bound_left_part)

    #TODO: сделать быстрее, а не перебор всех пар

    max_lcp = 0 # максимальная длина повтора
    coord1 = -1
    coord2 = -1

    for item in d.items():
        (left_indices_list, right_indices_list) = item[1]
        for index_left in left_indices_list:  # первая часть строки такая
            for index_right in right_indices_list: # сравниваем только со второй частью строки (после pos)

                two_parts_lcp = lcp(suff_matrix, index_left, index_right)
                c1 = beg_pos - lower_bound_left_part + index_left
                if end_pos is None:
                    c2 = beg_pos - lower_bound_left_part + index_right - 1
                else:
                    x = beg_pos -lower_bound_left_part + index_right - 1 - beg_pos
                    c2 = end_pos + x
                    # c2 = end_pos + beg_pos - lower_bound_left_part + index_right - 1
                if max_lcp < two_parts_lcp and two_parts_lcp < 51 and two_parts_lcp > 19 and c2 + two_parts_lcp - c1 < 5001:
                    max_lcp = two_parts_lcp
                    coord1 = c1
                    coord2 = c2



    # print(indexes_of_max_lcp, max_lcp)
    # print("Coordinates regarding pos = ", coord1, coord2)
    print("Direct repeat length =", max_lcp)
    print(coord1, '(' + str(beg_pos) + ', ' + str(end_pos) + ')', coord2)
    return (max_lcp, coord1, coord2)


# def find_direct_repeats(genome, beg_pos, end_pos):
#     joined_string, suff_matrix, lower_bound_left_part = create_suffix_matrix(genome, beg_pos, end_pos)
#     suff_array = suff_matrix[-1]
#
#     # index_list - суффиксы в лексикограяичесокм порядке
#     # pairs_array = массив пар, где в паре (массив с индексами слева от pos c 3меров в начале)
#     # справа массив с индексами спрва от pos с тем же 3мером вначале
#     index_list, pairs_array, d = create_index_massive_list(suff_array, joined_string, lower_bound_left_part)
#
#
#     max_lcp = 0 # максимальная длина повтора
#
#     # coord1 = -1
#     # coord2 = -1
#
#     # for item in d.items():
#     #     (left_indices_list, right_indices_list) = item[1]
#     #     for index_left in left_indices_list:  # первая часть строки такая
#     #         for index_right in right_indices_list: # сравниваем только со второй частью строки (после pos)
#     #
#     #             # two_parts_lcp = lcp(suff_matrix, index_left, index_right)
#     #             # c1 = beg_pos - lower_bound_left_part + index_left
#     #             #
#     #             # if end_pos is None:
#     #             #     c2 = beg_pos - lower_bound_left_part + index_right - 1
#     #             # else:
#     #             #     c2 = end_pos + (beg_pos - lower_bound_left_part) + index_right - 1
#     #             #
#     #             # #print(two_parts_lcp)
#     #             #
#     #             # if max_lcp < two_parts_lcp and c2 + two_parts_lcp - c1 < 5001 and two_parts_lcp < 51 and two_parts_lcp > 11:
#     #             #     max_lcp = two_parts_lcp
#     #             #     coord1 = c1
#     #             #     coord2 = c2
#     #
#     #             two_parts_lcp = lcp(suff_matrix, index_left, index_right)
#     #             if max_lcp < two_parts_lcp:
#     #                 max_lcp = two_parts_lcp
#     #                 indexes_of_max_lcp = (index_left, index_right)
#     #                 coord1 = index_left
#     #                 coord2 = index_right
#
#     max_lcp = 0  # максимальная длина повтора
#     indexes_of_max_lcp = (0, 0)
#
#     for item in d.items():
#         (left_indices_list, right_indices_list) = item[1]
#         for index_left in left_indices_list:  # первая часть строки такая
#             for index_right in right_indices_list:  # сравниваем только со второй частью строки (после pos)
#
#                 two_parts_lcp = lcp(suff_matrix, index_left, index_right)
#                 if max_lcp < two_parts_lcp:
#                     max_lcp = two_parts_lcp
#                     indexes_of_max_lcp = (index_left, index_right)
#
#     coord1 = beg_pos - lower_bound_left_part + indexes_of_max_lcp[0]
#     if end_pos is None:
#         coord2 = beg_pos - lower_bound_left_part + indexes_of_max_lcp[1] - 1
#     else:
#         coord2 = end_pos - lower_bound_left_part + indexes_of_max_lcp[1] - 1
#
#
#     # print(indexes_of_max_lcp, max_lcp)
#     # print("Coordinates regarding pos = ", coord1, coord2)
#     print("Direct repeat length =", max_lcp)
#     print(coord1, '(' + str(beg_pos) + ', ' + str(end_pos) + ')', coord2)
#     return (max_lcp, coord1, coord2)



def find_reverse_compl_repeats(genome, beg_pos, end_pos):
    joined_string, suff_matrix, lower_bound_left_part, upper_bound = create_reverse_compl_suffix_matrix(genome, beg_pos, end_pos)
    suff_array = suff_matrix[-1]

    # index_list - суффиксы в лексикограяичесокм порядке
    # pairs_array = массив пар, где в паре (массив с индексами слева от pos c 3меров в начале)
    # справа массив с индексами спрва от pos с тем же 3мером вначале
    index_list, pairs_array, d = create_index_massive_list(suff_array, joined_string, lower_bound_left_part)

    max_lcp = 0 # максимальная длина повтора

    coord1 = -1
    coord2 = -1

    for item in d.items():
        (left_indices_list, right_indices_list) = item[1]
        for index_left in left_indices_list:  # первая часть строки такая
            for index_right in right_indices_list: # сравниваем только со второй частью строки (после pos)

                two_parts_lcp = lcp(suff_matrix, index_left, index_right)

                c1 = beg_pos - lower_bound_left_part + index_left
                if end_pos is not None:
                    y = beg_pos - lower_bound_left_part + index_right -1 + two_parts_lcp - beg_pos
                    x = upper_bound - y
                    #print(y, upper_bound)
                    c2 = end_pos + x


                else:
                    x = beg_pos - lower_bound_left_part + index_right - 1 -beg_pos
                    y = upper_bound - x
                    c2 = beg_pos + y - two_parts_lcp

                if max_lcp < two_parts_lcp and c2 + two_parts_lcp - c1 < 5001 and two_parts_lcp < 51 and two_parts_lcp > 19:
                    max_lcp = two_parts_lcp
                    coord1 = c1
                    coord2 = c2

    print("Reverse repeat length =", max_lcp)
    print(coord1, '(' + str(beg_pos) + ', ' + str(end_pos) + ')', coord2)
    return (max_lcp, coord1, coord2)


def find_exact_repeats_near_pos(genome, beg_pos, end_pos=None):
    (max_lcp, coord1, coord2) = find_direct_repeats(genome, beg_pos, end_pos)
    (max_lcp_reverse, reverse_coord1, reverse_coord2) =find_reverse_compl_repeats(genome, beg_pos, end_pos)
    return max_lcp, max_lcp_reverse, coord1, coord2, reverse_coord1, reverse_coord2


def find_Elizavetas_repeat(genome_seq):
    repeat = "ACTCAAACTTCGCACACCTGAAACCCTTA"
    pos_founded_repeat_first = genome_seq.find(repeat)
    print("Pos1 Elizaveta's repeat =", pos_founded_repeat_first)
    # 1343644
    pos_founded_repeat_first = genome_seq[1343645:].find(repeat)
    print("Pos2 Elizaveta's repeat =", pos_founded_repeat_first + 1343644 + 1)
    # 1347291


    print(genome_seq[1345213: 1345213 + 14])
    print(genome_seq.find("TGCGAAGTTTGAGT"))
    print(genome_seq[1343645:].find("ACTCAAACTTCGCA") + 1343645 + 1)


def read_bed(filename):
    coordinates = []
    with open(filename, 'r') as f:
        for line in f:
            tmp = line.split()
            coord1 = int(tmp[1])
            coord2 = int(tmp[2])
            new_pair = (coord1, coord2)
            # убрали повторы, которые есть в бед файле
            if len(coordinates) == 0 or new_pair != coordinates[-1]:
                coordinates.append((coord1, coord2))

    return coordinates


def test():
    filename = "genome1.fasta"
    pos = 1346120
    genome_name, genome_seq = genome_read(filename)
    find_exact_repeats_near_pos(genome_seq, pos)
    find_Elizavetas_repeat(genome_seq)


def make_plots_length(statistic_array):

    d = {}

    for (length, dist) in statistic_array:
        try:
            d[length] += 1
        except:
            d[length] = 1


    max_length = max(list(d.keys()))

    x = [i for i in range(max_length + 1)]
    #x = sorted(list(d.keys()))
    #y = [d[i] for i in x]
    y_pos = np.arange(len(x))
    y = []
    for i in x:
       try:
           y.append(d[i])
       except:
           y.append(0)
    plt.bar(y_pos, y, align='center', alpha=0.5)
    plt.xticks(y_pos, x)
    plt.show()


def make_plots_dists(statistic_array):
    d = {}

    for (length, dist) in statistic_array:
        try:
            d[dist] += 1
        except:
            d[dist] = 1

    max_dist = max(list(d.keys()))

    x = sorted(list(d.keys()))
    y_pos = np.arange(len(x))
    y = [d[i] for i in x]

    plt.bar(y_pos, y, align='center', alpha=0.5)
    plt.xticks(y_pos, x)
    plt.title('DISTS')
    plt.show()



def run(filename, bedfile):


    genome_name, genome_seq = genome_read(filename)
    coords_list = read_bed(bedfile)
    print("Count of pos =", len(coords_list))

    w = open("repeats_" + bedfile, 'w')

    # массив, в котром в паре будут записаны пара длина повтора, расстояние от своей beg_pos
    statistic_array = []

    count_both_repeats = 0
    count_has_only_direct_repeats = 0
    count_has_only_reverse_repeats = 0
    count_without_repeats = 0
    for i, (beg_coord, end_coord) in enumerate(coords_list):
        print(i)
        #TODO птоврты только больше длины 3
        direct_repeat_length, reverse_repeat_length, drcoor1, drcoor2,\
            rrcoord1, rrcoord2 = find_exact_repeats_near_pos(genome_seq, beg_coord, end_coord)

        print("OUR= ",drcoor1, drcoor2, "REV",  rrcoord1, rrcoord2)

        # TODO точно ли минимумы
        statistic_array.append((direct_repeat_length, min(beg_coord - drcoor1, drcoor2-beg_coord)))
        statistic_array.append((reverse_repeat_length, min(beg_coord-rrcoord1, rrcoord2 - beg_coord)))

        if direct_repeat_length != 0:
            w.write("Direct" +  ' ' + str(drcoor1) + ' ' + str(drcoor2) + ' ' + str(direct_repeat_length) +'\n')
        if reverse_repeat_length != 0:
            w.write("Reverse" + ' ' + str(rrcoord1) + ' ' + str(rrcoord2) + ' ' + str(reverse_repeat_length) + '\n')




        if direct_repeat_length > 3 and reverse_repeat_length > 3:
            count_both_repeats += 1
        elif direct_repeat_length > 0:
            count_has_only_direct_repeats += 1
        elif reverse_repeat_length > 0:
            count_has_only_reverse_repeats += 1
        else:
            count_without_repeats += 1

    w.close()
    print("Count of pos =", len(coords_list))
    print("Count both repeats =", count_both_repeats)
    print("Count only direct =", count_has_only_direct_repeats)
    print("Count reverse direct =", count_has_only_reverse_repeats)
    print('Count without repeats =', count_without_repeats)


    # make_plots_dists(statistic_array)


if __name__ == "__main__":
    genome = sys.argv[1]
    bedfile = sys.argv[2]
    run(genome, bedfile)

