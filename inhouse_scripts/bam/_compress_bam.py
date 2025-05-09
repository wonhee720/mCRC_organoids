import re
import pysam 


CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
CIGAR2CODE = dict([y, x] for x, y in enumerate("MIDNSHP=XB"))

def reads_to_orientation(v):
    if v[0].is_reverse: 
        v_0 = 'R'
    else: 
        v_0 = 'F'
    if v[1].is_reverse: 
        v_1 = 'R'
    else: 
        v_1 = 'F'
    if v[0].is_read1:
        v_0 += '1'
    elif v[0].is_read2:
        v_0 += '2'
    if v[1].is_read1:
        v_1 += '1'
    elif v[1].is_read2:
        v_1 += '2'
    return ''.join(sorted([v_0, v_1]))


def bam_to_reads(input_bam, chromosome=None):
    reads = defaultdict(list)
    for read in input_bam.fetch(contig=chromosome):
        reads[read.query_name].append(read)
    return reads


def seq_to_matrix(sequence, length=151):
    base_to_number = {'A':0, 'C':1, 'G':2, 'T':3}
    matrix = np.zeros((4, length))
    for j, i in enumerate(sequence):
        try:
            matrix[base_to_number[i], j] = 1
        except:
            continue
    return matrix


def matrix_to_seq(matrix):
    number_to_base = {0:'A', 1:'C', 2:'G', 3:'T'}
    seq = ""
    for j in range(matrix.shape[1]):
        try:
            seq += number_to_base[np.argmax(matrix[:,j])]
        except:
            continue
    return seq
            

def cigartuples_to_number(cigar_tuples, start=True):
    if cigar_tuples is None:
        return 0
    cigar_list = list(zip(*cigar_tuples))
    index = cigar_list[0].index(0)
    if start:
        return sum(cigar_list[1][:index])
    else:
        return sum(cigar_list[1][index:])


def collapse_umis_to_bam(k, v, read_orientation, output_bam):
    number_of_reads = len(v[read_orientation])
    if number_of_reads:
        read_name = ':'.join([str(x) for x in k]) + ':'+read_orientation+':' + str(number_of_reads)
        first_query_length = max([x[0].query_length for x in v[read_orientation]]) # v['F1R2'][0][0].query_length
        second_query_length = max([x[1].query_length for x in v[read_orientation]]) # v['F1R2'][0][1].query_length
        r1 = copy.deepcopy(v[read_orientation][0][0])
        r2 = copy.deepcopy(v[read_orientation][0][1])
        r1_sequence_mtx = np.zeros((4, first_query_length))
        r2_sequence_mtx = np.zeros((4, second_query_length))
        
        for reads_f1r2 in v[read_orientation]:
            r1_sequence = reads_f1r2[0].query_sequence
            r2_sequence = reads_f1r2[1].query_sequence
            r1_qualities = np.zeros(first_query_length)
            r2_qualities = np.zeros(second_query_length)
            r1_qualities[0:len(reads_f1r2[0].query_qualities)] = reads_f1r2[0].query_qualities
            r2_qualities[0:len(reads_f1r2[1].query_qualities)] = reads_f1r2[1].query_qualities
            r1_sequence_mtx += seq_to_matrix(r1_sequence, length=first_query_length) * r1_qualities
            r2_sequence_mtx += seq_to_matrix(r2_sequence, length=second_query_length) * r2_qualities
            
        r1_qualities = array.array('B', (np.max(r1_sequence_mtx, axis=0) / number_of_reads).astype(int))
        r2_qualities = array.array('B', (np.max(r2_sequence_mtx, axis=0) / number_of_reads).astype(int))
        r1_sequence = matrix_to_seq(np.where(r1_sequence_mtx == r1_sequence_mtx.max(axis=0, keepdims=True), 1, 0).astype(int))
        r2_sequence = matrix_to_seq(np.where(r2_sequence_mtx == r2_sequence_mtx.max(axis=0, keepdims=True), 1, 0).astype(int))

        r1.query_name = read_name
        r2.query_name = read_name
        r1.query_sequence = r1_sequence
        r2.query_sequence = r2_sequence
        r1.query_qualities = r1_qualities
        r2.query_qualities = r2_qualities
        r1.cigarstring = None
        r2.cigarstring = None

        v_collapse = [r1, r2]
        for v in v_collapse:
            output_bam.write(v)
        return
    else:
        return
            