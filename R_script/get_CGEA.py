import pickle
import os
import CGEA
import pandas as pd
from CGEA.lib.statistics import generateBackgroundModel, scoreRankedLists, generatePVals
from  CGEA.lib.featureListProcessing import formatGeneLists
from time import time
from multiprocessing import Pool

def process_res(result):                     
    res = []
    for item in result['Compound Connectivity Scores']:
           res.append((item['Name'], item['DrugBank'], item['Adjusted P-Value'], item['P-Value'], item['Score']))
    return res
    
def featureBasedCGEA_simpliefied(upgenes, downgenes, cmap, numWorkers = 6, numPermutations = 1000, pool = None, alpha = 0.05):
    geneLists = formatGeneLists([upgenes, downgenes])
    
    [backgroundModel, backgroundScores, scaler] = generateBackgroundModel(geneLists, cmap, 
                                                                          numWorkers = numWorkers, numPermutations = numPermutations)
    [drugIndex, scores, index] = scoreRankedLists(geneLists, cmap, numWorkers = numWorkers, pool = pool)
    
    [pvals, adjpvals, scores] = generatePVals(backgroundModel, scores, backgroundScores = backgroundScores, scaler = scaler, 
                                              alpha = alpha, numWorkers = numWorkers, pool = pool)

    drugAnnotations = {}
    for drug in drugIndex.keys():
        index = drugIndex[drug]
        drugAnnotations[drug]={'Score':scores[index],
                               'P-Value':pvals[index],
                               'Adjusted P-Value':adjpvals[index]
                               }

    return drugAnnotations

def par_proc(var_in):
    tissue, key, gs, cmap = var_in[0], var_in[1], var_in[2],  var_in[3]
    res = []
    try:
        res = featureBasedCGEA_simpliefied(gs, [], cmap, numWorkers = 1, numPermutations = 1000, pool = None, alpha = 0.05)
    except:
        print tissue, key
    return key, res
    
data = CGEA.lib.CGEAdata.CGEAdata()
cmap = data.get('cmap')

path_wgcna = './WGCNA/'
path_cgea = './CGEA/'
pool = Pool(20)

if __name__ == '__main__':
    fns = sorted(os.listdir(path_wgcna))
    
    for tissue in fns:
        if os.path.exists(path_cgea + ('%s.p' % tissue)):
            continue
        
        try:
            df_wgcna = pd.read_csv(path_wgcna + '%s/gene_module.csv' % tissue, header = None)
        except:
            continue
            
        with open('../data/gene_info/mart_export.txt', 'rb') as f:
            data = f.readlines()

        data = [item.strip().split(',') for item in data[1:]]

        entrez2ensembl = {item[2]:item[0] for item in data if item[2] != ''}
        ensembl2entrez = {g2:g1 for g1, g2 in entrez2ensembl.items()}

        mod2genes ={k:[ensembl2entrez[g] for g in item[0] if g in ensembl2entrez] for k, item in df_wgcna.groupby(1) if k != 'grey'}

        inputs = [(tissue, key, gs, cmap) for key, gs in mod2genes.items()]
        res = pool.map(par_proc, inputs)
        cgea_res = {k:r for k, r in res}

        with open(path_cgea + ('%s.p' % tissue), 'wb') as f:
            pickle.dump(cgea_res, f)