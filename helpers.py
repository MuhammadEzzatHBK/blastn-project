"""
Created on Wed Dec 16 15:21:49 2020

@author: Sara Ahmad Maher
         Muhammad Ayman Ezzat
         Salma Muhammad Saeed
         Salma Muhammad Kamel
         Youmna Magdy Abdullah
         Sahar Saber Ibrahim
"""
import numpy as np
""" Extracts every possible kmer from a sequence.

    Args:
        sequence(str): Sequence to extracts kmers from
        k(int): Word length
              
    Returns:
        kmers(list)(str):Every possible sequence kmer
 """
def extract_kmers(sequence,k):
    kmers=[]
    for i in range(len(sequence)-k+1):
        kmers.append([sequence[i],sequence[i+1],sequence[i+2]])
    return kmers


""" Encodes a certain kmer into a key.

    Args:
        kmer(str): kmer to encode
              
    Returns:
        code(int): kmer key
        
    Note: May not be needed.
 """
def encode_kmers(kmer):
    encoder = {'A':0,'T':3,'G':2,'C':1}
    code = 0
    for i in range(len(kmer)):
        code += encoder[kmer[i]] * 4**(len(kmer)-(i+1))
    return code

""" Preforms naive ungapped alignment between two kmers.

    Args:
        kmer1(str): First Sequence
        kmer2(str): Second Sequence
        match_score(float): Score added on match occurunce in ungapped alignment
        mismatch_score(float) : Score added on mismatch occurunce in ungapped alignment
              
    Returns:
        (float) Alignment Score
 """
def ungapped_alignment(kmer1,kmer2,match_score,mismatch_score):
    matrix = np.zeros((len(kmer1)+1,len(kmer2)+1))
    for i in range(1,len(kmer1)+1):
        for j in range(1,len(kmer2)+1):
            if(kmer1[i-1]==kmer2[j-1]):
                matrix[i][j]= matrix[i-1][j-1] + match_score
            else:
                matrix[i][j]= matrix[i-1][j-1] + mismatch_score
    return matrix[len(kmer1)][len(kmer2)]

