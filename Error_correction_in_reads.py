from Bio import SeqIO

HmDistance = 1
complement_table = {"A":"T","T":"A","G":"C","C":"G"}
file = "F:\\avans\\stage erasmus\\rosalin\\rosalind_corr.txt"


def reverse_compliment(seq):

    ''''
    calculates the reverse compliment from the given sequention

    parameters:
    ----------
    seq: str
        sequention of DNA

    :return
    --------
    reverse compliment of given DNA sequention
    '''

    rev_complement = ""
    for N in (seq[::-1]):
        rev_complement = rev_complement + complement_table[N]
    return(rev_complement)

def read_in_fasta(file):
    ''''
    reads in the given fasta file

    :param:
    file: str
        path to fasta file

    :return:
        sequentions:list
            list of sequentions read in of the fasta file
    '''
    sequentions = []
    with open(file, 'r') as fa:
        for seq_record in SeqIO.parse(fa, 'fasta'):
            sequentions.append(str(seq_record.seq))
    return(sequentions)

#retrieve the sequentions
sequentions = read_in_fasta(file)
correct = []

#filtering the correct sequentions
for i in range(len(sequentions)):
    for ii in range(len(sequentions)):
        #little check that they dont select it self in the list (based on numbers)
        if i != ii:
            if sequentions[i] == sequentions[ii] or sequentions[i] == reverse_compliment(sequentions[ii]):
                correct.append(sequentions[i])

#select mutated
errors = []
for i in sequentions:
    if i not in correct:
        errors.append(i)

#remove duplicated correct sequentions to reduce calc time.
correct = list(dict.fromkeys(correct))


corrected = []
#loops through all sequentions with an error
for i in errors:
    #check every correct sequention
    for ii in correct:

                #check on normal
                desimilar = 0
                for n in range(len(i)):
                    if i[n] != ii[n]:
                        desimilar += 1
                if desimilar == HmDistance:
                    corrected.append('{0}->{1}'.format(i,ii))
                    break

                #check on reverse compliment
                rev_ii = reverse_compliment(ii)
                desimilar = 0
                for n in range(len(i)):
                    if i[n] != rev_ii[n]:
                        desimilar += 1

                if desimilar == HmDistance:
                    corrected.append('{0}->{1}'.format(i, rev_ii))
                    break

for i in corrected:
    print(i)
