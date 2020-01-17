from multiprocessing import Pool
import pandas as pd
import numpy as np
from scipy.stats import linregress
import os
from itertools import product
import gzip
import time
import argparse

def normalize_1d(x):
    ## x: (n, ) np array
    return (x - np.mean(x)) / (np.std(x) + 1e-5)
        
def par_proc(var_in):
    file_in, file_out, dosage_file, refIID = var_in
    
    tic = time.clock()
    
    dosage_dict, rsid2allel = read_dosage_files(dosage_file)
    
    ## read trait
    with open(file_in) as f:
        tmp = [item.strip().split('\t') for item in f]
        traitIID = [item[1] for item in tmp[1:]]
        trait = [float(item[2]) for item in tmp[1:]]

    tmp = set(traitIID)
    
    ## find common IIDs
    id2IID = [d for d in refIID if d in tmp]
    IID2id = {d:i for i, d in enumerate(id2IID)}
    
    ## trait vector
    trait_vec = [0.0] * len(id2IID)
    for d, t in zip(traitIID, trait):
        if d in IID2id:
            trait_vec[IID2id[d]] = t

    ## indeice for dosage vector 
    refIID2id = {d:i for i, d in enumerate(refIID)}
    ind_ref = np.array([refIID2id[d] for d in id2IID])

    
    fs = []
    for rsid, row in dosage_dict.items():
        dosage_vec = row[ind_ref].copy()
        MAF = dosage_vec.mean() / 2.0 
        if MAF < 0.01 or MAF > 0.99:
            continue

        dosage_vec = normalize_1d(dosage_vec)
        
        ## linear regression
        slope, intercept, r_value, p_value, std_err = linregress(dosage_vec, trait_vec)
        a1, a2 = rsid2allel[rsid]
        fs.append((rsid, a1, a2, str(slope), str(p_value)))
    
    ## output files
    with gzip.open(os.path.join(file_out), 'wb') as writer:
        writer.write('\t'.join(['SNP', 'A1', 'A2', 'BETA', 'P']) + '\n')

        for item in fs:
            writer.write('\t'.join(item) + '\n')

    toc = time.clock()
    print toc - tic
            
def read_dosage_files(file_in):
    rsid2allel = {}
    dosage_dict = {}
    with gzip.open(file_in, 'rb') as f:
        for i, item in enumerate(f):
            tmp = item.split()
            dos_vec = map(float, tmp[6:])
            rsid = tmp[1]
            rsid2allel[rsid] = (tmp[3], tmp[4])
            dosage_dict[rsid] = np.array(dos_vec)
            #if i > 100:
            #    break
                
    return dosage_dict, rsid2allel

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser = argparse.ArgumentParser(description='parser')
    parser.add_argument('--folder_dosage', action="store")
    parser.add_argument('--folder_experiment', action="store")
    parser.add_argument('--n_cores', action="store", default=25)
    args = parser.parse_args()

    
    folder_to_dosages = args.folder_dosage ## assumes there are files: samples.txt and chr(1-22)_1kg_phase3_dosages.tsv.gz
    folder_to_exp = args.folder_experiment
    p = Pool(int(args.n_cores))
    

    ## implied path for traits
    folder_to_traits = os.path.join(folder_to_exp, 'traits')
    
    ## implied path for snp_associations, will create if not exist
    folder_to_snp_ass = os.path.join(folder_to_exp, 'snp_associations')
    if not os.path.exists(folder_to_snp_ass):
        os.makedirs(folder_to_snp_ass)   
        
    ## read IID
    refIID_files = folder_to_dosages + 'samples.txt'
    with open(refIID_files) as f:
        refIID = [item.strip().split(' ')[1] for item in f]    
    
    dosage_files = [folder_to_dosages  + 'chr%d_1kg_phase3_dosages.tsv.gz' % i for i in range(1, 23)]

    for _, _, files in os.walk(folder_to_traits, topdown=False):
        files = sorted(files)
        
        for name in files:
            
            file_in = os.path.join(folder_to_traits, name)
            folder_out = os.path.join(folder_to_snp_ass, name.split('.')[0])
            
            if not os.path.exists(folder_out):
                os.makedirs(folder_out)
            
            vars_in = []
            for i in range(22):
                file_out = os.path.join(folder_out, 'chr%i.tsv.gz' % (i + 1))
                vars_in.append((file_in, file_out, dosage_files[i], refIID))
                #vars_in.append((file_in, file_out, dosage_dict, rsid2allel, refIID ))

            p.map(par_proc, vars_in)
