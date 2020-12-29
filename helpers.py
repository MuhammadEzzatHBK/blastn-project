"""
Created on Wed Dec 16 15:21:49 2020

@author: Sara Ahmed Maher
         Salma Muhammad Saied
         Salma Muhammad Kamel
         Youmna Magdy Abdullah
         Sahar Saber Ibrahim     
         Muhammad Ayman Ezzat

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
    kmer = ""
    for i in range(len(sequence)-k+1):
        kmer = sequence[i]
        for j in range(1,k):
            kmer+=sequence[i+j]      
        kmers.append(kmer)
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
    scores = np.zeros(len(kmer1))
    for i in range(len(kmer1)):
        if(kmer1[i-1]==kmer2[i-1]):
            scores[i] = scores[i-1]+match_score
        else:
            scores[i] = scores[i-1]+mismatch_score
    return scores[-1]

""" The core process of extention through Smith-Waterman algorithm.
     
    Args:
        hssp_table(data.frame) : A dataframe of hssps, the return of find_hssp function.
        querry_seq(str) : The querry sequence
        db_seq(str) : The database sequence
        k(int): Word length
        match_score(int): Socre given on match occurunce
        mismatch_score(int) : Socre given on mismatch occurunce
        gap_score(int) : Socre given on gap placment
        extention_threshold(int) : The threshold alignment score that if an extention falls
                                   behind it, the extention stops
    Returns : 
        raw_score(int): Raw alignment score between two extended kmers
        querry_alignment : querry kmer after extention
        db_alignment : db kmer after extention
"""
def SWM(querry_seq,db_seq,querry_index,db_index,k,
        match_score,mismatch_score,gap_score,
        extend_without_checking, extention_threshold,score):
   querry_alignment = ""
   db_alignment = ""
   raw_score = 0
   matrix = np.zeros((len(querry_seq),len(db_seq)))
   matrix[querry_index+k-1][db_index+k-1] = score
   resume = True
   row = len(querry_seq) - 1
   col = len(db_seq) - 1
   #Start building table from the kmer position
   start_checking =querry_index+k+extend_without_checking
   for i in range(querry_index+k,len(querry_seq)):
       if(not resume):
           break
       for j in range(db_index+k,len(db_seq)):
           R = 0
           C = 0
           if(querry_index+k+1-i >= len(querry_seq)):   
               break
           if(querry_seq[i] == db_seq[j]):
               D = matrix[i-1][j-1] + match_score
           else:
               D = matrix[i-1][j-1] + mismatch_score
               R = matrix[i][j-1] + gap_score
               C = matrix[i-1][j] + gap_score
               
           matrix[i][j] = np.max([D,R,C])
           raw_score = matrix[i][j] 
          #Break if we drop below threshold
           if(i>= start_checking):
               if(raw_score < extention_threshold):
                   resume = False
                   row = i
                   col = j
                   break
               else:
                   extention_threshold = raw_score 
    
   
   while(row > querry_index & col>db_index):
        if((matrix[row][col] == matrix[row-1][col-1]+match_score) | (matrix[row][col] == matrix[row-1][col-1]+mismatch_score) ):
            querry_alignment += querry_seq[row-1]
            db_alignment += db_seq[col-1]
            row-=1
            col-=1
        elif(matrix[row][col] == matrix[row-1][col]):
            querry_alignment += querry_seq[row]
            db_alignment += "_"
            row -=1
        else:
            db_alignment += db_seq[col]
            querry_alignment += "_"
            col -=1            
           
   querry_alignment = querry_alignment[::-1]
   db_alignment = db_alignment[::-1]                
   querry_alignment = querry_seq[querry_index:querry_index+k+1] + querry_alignment
   db_alignment = db_seq[db_index:db_index+k+1] + db_alignment
    
   return raw_score,querry_alignment,db_alignment
                 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
