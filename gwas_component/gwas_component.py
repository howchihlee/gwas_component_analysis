import os
import pandas as pd
import numpy as np
from collections import Counter
import pickle
import gzip
from scipy.stats import norm
from statsmodels.stats import multitest as multi

def metaxcan_file_to_gene2z(fn):
    df = pd.read_csv(fn)
    gene2z = {}
    for g, z in zip(df.gene, df.zscore):
        gene2z[g.split('.')[0]] = z
    return gene2z

def compute_gwas_z(pca_genes, pca_weight, pca_sig, gene2zscore, gene2var):
    ## inputs:
    ## df_spca: a (lantent factor, gene) numpy array of the lantent factor x gene matrix
    ## sig_gtex: a (gene, ) numpy array of the per-gene variance from GTEX
    ## df_cov: a (gene, gene) numpy array of the gene-gene covariance
    ## z_g: a (gene, ) numpy array of per-gene GWAS associations from s-predixscan
    ## var_g: a (gene, ) numpy array of per-gene GWAS variance from s-predixscan

    z_lantent = 0.0
    for i, g in enumerate(pca_genes):
        if g in gene2zscore and g in gene2var:
            w = pca_weight[i]
            z = gene2zscore[g]
            s = np.sqrt(gene2var[g])

            z_lantent += w*z*s
    z_lantent = z_lantent / pca_sig
    return z_lantent
    
## constants:
## path to model 
## list of tissues
_DIR_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__)))                 
_MODEL_PATH = os.path.abspath(os.path.join(_DIR_PATH, './eigengene_model/'))
_GENEVAR_PATH = os.path.abspath(os.path.join(_DIR_PATH, './gene_info/gene_var/'))
_GENE_PATH = os.path.abspath(os.path.join(_DIR_PATH, './gene_info/model_genes/'))
_ANNOTATION_PATH = os.path.abspath(os.path.join(_DIR_PATH, './module_annotation/'))

def _get_model_name():
    return [s.split('.')[0] for s in os.listdir(_MODEL_PATH)]

def _get_modelgene(tissue0):
    with gzip.open(os.path.join(_GENE_PATH, '%s-model_genes.txt.gz' % tissue0), 'rb') as reader:
        id2gene = [g.strip() for g in reader.readlines()]
    return id2gene

def _get_gene2var(tissue0):
    with gzip.open(os.path.join(_GENEVAR_PATH, '%s-gene2var.csv.gz' % tissue0), 'rb') as reader:
        data = reader.readlines()
        gene2var = {}
        for g in data:
            g, v = g.strip().split(',')
            gene2var[g.split('.')[0]] = float(v)
    return gene2var

def _get_eigengene_models(tissue0):
    with open(os.path.join(_MODEL_PATH, '%s.p' % tissue0), 'rb') as reader:
        pca_models = pickle.load(reader)
    return pca_models    
    
def _get_annotation(tissue0):
    try:
        with open(os.path.join(_ANNOTATION_PATH, '%s.p' % tissue0), 'rb') as reader:
            annotations = pickle.load(reader)  
        return annotations
    except: 
        return {}
    
class Coexpression_Model():
    def __init__(self, tissue0):
        self.tissue = tissue0
        self.id2gene = _get_modelgene(tissue0)
        self.gene2var = _get_gene2var(tissue0)
        pca_models = _get_eigengene_models(tissue0)
        self.pca_models = [(k, pca_genes, pca_weight, sig_1k) for k, pca_genes, pca_weight, sig_1k in pca_models if k != 'grey']
        self.annotations = _get_annotation(tissue0)
        self.is_computed = False
        self.res_z = None
        self.res_p = None
        self.res_adjp = None
        self.adj_asso_components = None
    
    def _print_not_fitted_message(self):
        print('plz run compute_z first.')
        return
    
    def compute_z(self, gene2zscore, fdr = 0.05):
        res_z = []
        res_p = []
        
        for k, pca_genes, pca_weight, sig_1k in self.pca_models:
            z = compute_gwas_z(pca_genes, pca_weight, sig_1k, gene2zscore, self.gene2var)
            res_z.append(z)
            res_p.append(norm.sf(abs(z)) * 2)
        
        res_z = np.array(res_z)
        res_p = np.array(res_p)
        _, res_adjp, _, _ = multi.multipletests(res_p)
        
        self.res_z = res_z
        self.res_p = res_p
        self.res_adjp = res_adjp
        self.adj_asso_components = np.where(res_adjp < fdr)[0]
        self.is_computed = True
        
        return res_z, res_p, res_adjp
        
    def get_driver_genes(self):
        if not self.is_computed:
            self._print_not_fitted_message()
            return 
        
        res_all = []
        for i in self.adj_asso_components:
            k = self.pca_models[i][0]
            if k not in self.annotations:
                continue
                
            if 'driver_genes' in self.annotations[k]:
                res_all += [(self.tissue, g) for g in self.annotations[k]['driver_genes']]
        return res_all    
    
    def get_associated_components(self):
        if not self.is_computed:
            self._print_not_fitted_message()
            return 
        
        res_all = []
        for i in self.adj_asso_components:
            k, pca_genes, _, _ = self.pca_models[i]
            res_all += [(k, pca_genes)]
        return res_all   
        