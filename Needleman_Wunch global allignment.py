import numpy as np

##############
# parameters #
##############

seqA = "GATTACA"
seqB = 'GTCGACGCA'

match = 1
mismatch = -1
gap = -2

def generate_matrix(match,mismatch,gap,seqA,seqB):

    '''
    Generates the matrix for the Needleman Wunch algorithm
    - creates a matrix of matches
    - create a basic gap matrix

    -combines each other to recalculate the correct main matrix

    :param match: int
        reward for a match
    :param mismatch: int
        penalty for mismatch
    :param gap: int
        penalty for a gap
    :param seqA: str
        Sequence to be alligned
    :param seqB: str
        sequence to be alligned
    :return:
        match_check : df
            dataframe containing whether a sequence is a match or not
        main_matrix: df
            dataframe containing the calculated matrix for the Needleman_wunsch algorithm

    '''

    #Creat empty Matrix
    main_matrix = np.zeros((len(seqA)+1,len(seqB)+1))
    match_check = np.zeros((len(seqA),len(seqB)))
    
    #fill matrix with match and mismatch
    for N in range(len(seqA)):
        for NN in range(len(seqB)):
            if seqA[N] == seqB[NN]:
                match_check[N][NN] = match
            else:
                match_check[N][NN] = mismatch
    
    #create in the mainmatrix the gap value (upper row & left column)
    for i in range(len(seqA)+1):
        main_matrix[i][0] = i*gap
    for i in range(len(seqB)+1):
        main_matrix[0][i] = i*gap
    
    #fill the main matrix with the mismatch
    for A in range(1,len(seqA)+1):
        for B in range(1,len(seqB)+1):
            #value will be assigned which one has the highest value
            main_matrix[A][B] = max(main_matrix[A-1][B-1]+match_check[A-1][B-1], #check diagonal left
                                    main_matrix[A-1][B]+gap,                     #check left
                                    main_matrix[A][B-1]+gap)                     #check up

    return(main_matrix,match_check)

def global_allignment(main_matrix,match_check,seqA,seqB):

    '''
    starts generating sequence allignment based on the main matrix values
    process stops if both sequences have been processed
    -starts in the corner Right below and goes to left upper corner (not the gap calc)

    -calculates the best way back based on if the value from a location is the same as previous
        with correction of the gap panalty or mismatch,

        -incase it isnt will check if can go to left value is calculates as tested value

        -otherwise will go up


    :param main_matrix: df
            dataframe containing the calculated matrix for the Needleman_wunsch algorithm
    :param match_check: df
            dataframe containing whether a sequence is a match or not
    :param seqA: str
        given sequention
    :param seqB: str
        given sequention
    :return:
        allignt sequention : str
            sequence allignt to the other sequence
    '''

    aligned_1 = ""
    aligned_2 = ""
    IA= len(seqA)
    IB = len(seqB)

    print(main_matrix)
    print(match_check)
    while (IA > 0 and IB > 0):


        if (IA > 0 and IB > 0 and main_matrix[IA][IB] == main_matrix[IA - 1][IB - 1] + match_check[IA - 1][
            IB - 1]):

            aligned_1 = seqA[IA - 1] + aligned_1
            aligned_2 = seqB[IB - 1] + aligned_2

            IA -= 1
            IB -= 1

        elif (IA > 0 and main_matrix[IA][IB] == main_matrix[IA - 1][IB] + gap):
            aligned_1 = seqA[IA - 1] + aligned_1
            aligned_2 = "-" + aligned_2

            IA -=  1
        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = seqB[IB - 1] + aligned_2

            IB -= 1
    return(aligned_1, aligned_2)

main_matrix,match_check = generate_matrix(match,mismatch,gap,seqA,seqB)
aligned_1, aligned_2 = global_allignment(main_matrix,match_check,seqA,seqB)

#result
print(aligned_1)
print(aligned_2)
