"""
Created on Wed Dec 16 15:21:49 2020

@author: Sara Ahmad Maher
         Muhammad Ayman Ezzat
         Salma Muhammad Saeed
         Salma Muhammad Kamel
         Youmna Magdy Abdullah
         Sahar Saber Ibrahim
"""
import pandas as pd
import math
runfile('helpers.py')

""" Finds the high scoring segment pairs.

    Args:
        querry_seq(str) : The querry sequence we search for
        db_seq(str) : One of the source database's sequences
        k(int) : Word length
        match_score(float) : Score added on match occurunce in ungapped alignment
        mismatch_score(float) : Score added on mismatch occurunce in ungapped alignment
        hssp_threshold(float) : The score threshold at which we decide whether the 
                         considered kmer is a hssp or not
       
    Returns:
        hssp_table(data.frame) : The high scoring segment pairs table, contains the kmer,
                     the starting index of querry kmer, the starting index of
                     db_kmer & the ungapped alignment score.
"""
def find_hssps(querry_seq,db_seq,k,match_score,mismatch_score,hssp_threshold):   
    querry_kmers = extract_kmers(querry_seq,k)
    db_kmers = extract_kmers(db_seq,k)
    kmers =[]
    q_indicies=[]
    db_indicies=[]
    scores=[]
    for i in range(len(querry_kmers)):
        for j in range(len(db_kmers)):
            score = ungapped_alignment(querry_kmers[i],db_kmers[j],match_score,mismatch_score)
            if(score>=hssp_threshold):
                kmers.append(db_kmers[j])
                q_indicies.append(i)
                db_indicies.append(j)
                scores.append(score)
    data = {'db_kmer':kmers,
            'querry_index':q_indicies,
            'db_index':db_indicies,
            'score': scores}
    hssp_table = pd.DataFrame(data,columns = ["db_kmer","querry_index","db_index","score"])        
    return hssp_table



def seed_extend():
    return 0



""" Calculates the expected value E.
    Args:
        k(float) : 
        m(int) : db_length
        n(int) : querry sequence length
        lamb(float) : 
        raw_score(float) : kmer raw score after extention
        
    Returns:
        E(float)
"""
def final_score(k,m,n,lamb,raw_score):
    return k*m*n*math.exp(-lamb*raw_score)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    