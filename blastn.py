"""
Created on Wed Dec 16 15:21:49 2020

@author: Sara Ahmed Maher
         Salma Muhammad Saied
         Salma Muhammad Kamel
         Youmna Magdy Abdullah
         Sahar Saber Ibrahim     
         Muhammad Ayman Ezzat

"""
exec(open('main_processes.py').read())
def blastn(filepath,querry_seq,k,
           match_score,mismatch_score,gap_score,
           hssp_threshold,extend_without_checking,extention_threshold):
    db_sequences = import_data(filepath)
    for i in range(len(db_sequences)):
        print("Sequence number",i+1,"results : "+ db_sequences[i] +"\n")
        print("============================================================\n")
        print(blast_pipeline(querry_seq,db_sequences[i],k,match_score,mismatch_score,gap_score,hssp_threshold,extend_without_checking,extention_threshold))
        print("============================================================\n")

blastn('TestCases/db.db','ACAATTC',3,2,-1,-1,2,1,4)