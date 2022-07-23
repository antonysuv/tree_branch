#!/usr/bin/env python3
import pandas
from itertools import product
import sys, argparse, os
import numpy as np
from math import log, ceil
from scipy.stats import multinomial, chi2
from math import factorial

###Convert TRAIN TEST to numeric and save as nupy array
def get_n_taxa(aln_file):
    aln_file_open=open(aln_file)
    all_taxa = [ str(r) for r in aln_file_open if r[0]==">"]
    n_taxa = len(set(all_taxa))
    return(n_taxa)


#Read FASTA convert to numeric
def fasta_pars(aln_file,seq_number,Lmax):
    aln=open(aln_file)
    dic={'A':'0','T':'1','C':'2','G':'3','-':'4'}
    matrix_out=[]
    fasta_dic={}
    for line in aln:
        if line[0]==">":
            header=line[1:].rstrip('\n').strip()
            fasta_dic[header]=[]
        elif line[0].isalpha() or line[0]=='-':
            for base, num in dic.items():
                line=line[:].rstrip('\n').strip().replace(base,num)
            line=list(line)
            line=[int(n) for n in line]
            #Mkae all matrices of equal length for CNN +[-15]*(Lmax-len(line)) 
            fasta_dic[header]+=line+[-15]*(Lmax-len(line)) 
            if len(fasta_dic)==seq_number:
                taxa_block=[]
                for taxa in sorted(list(fasta_dic.keys())):
                    taxa_block.append(fasta_dic[taxa.strip()])
                fasta_dic={}
                matrix_out.append(taxa_block)
    return(np.array(matrix_out))

#Read TRAIN and TEST datasets to equalize sizes 
def tv_parse(train,test):
    seq_number = get_n_taxa(test)
    print("N taxa: "+str(seq_number))
    tr=open(train)
    te=open(test)
    LT=max([len(r.strip()) for r in tr])
    print("TRAIN largest alignment: "+str(LT))
    LTE=max([len(r.strip()) for r in te])
    print("TEST largest alignment: "+str(LTE))
    Lmax=max([LT]+[LTE])
    tr=fasta_pars(train,seq_number,Lmax)
    te=fasta_pars(test,seq_number,Lmax)
    print("N TRAIN alignments: "+str(len(tr)))
    print("N TEST alignments: "+str(len(te)))
    return tr, te   
    
def main():
    parser = argparse.ArgumentParser(description='fasta2numeric conversion Ready for Keras')
    parser.add_argument( '--tr', help = "Train dataset in FASTA",dest='TRAIN')
    parser.add_argument( '--te', help = "Test dataset in FASTA",dest='TEST')
    args = parser.parse_args()
    
    print("Reading input")
    train_data1, test_data1 = tv_parse(args.TRAIN,args.TEST)    
    #Reshape for Keras
    train_data1=train_data1.reshape(train_data1.shape[0],train_data1.shape[1],train_data1.shape[2],1)
    test_data1=test_data1.reshape(test_data1.shape[0],test_data1.shape[1],test_data1.shape[2],1)
    np.save("TRAIN", train_data1)
    np.save("TEST", test_data1)
     
if __name__ == "__main__":
    main()