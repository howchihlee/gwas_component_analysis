import pandas as pd
import os
from multiprocessing import Pool
        
def run_metaxcan(cmd):
    return os.system(cmd)

if __name__ == '__main__':
    p = Pool(20)
    cmds = []

    gwas_path = os.path.abspath('./single_component/snp_associations/')
    files = [f for f in os.listdir(gwas_path)]
    metaxcan_file = os.path.abspath('../metaxcan/MetaXcan.py') 

    for tissue in ['Whole_Blood']:
        for exp in files:

            fn0 = os.path.join(gwas_path, exp)
            fn1 = os.path.abspath('./single_component/metaxcan/%s_%s.csv' % (exp, tissue))
            cmd = 'python %s ' % metaxcan_file
            cmd += '--model_db_path %s ' % ('/data2/users/hclee/data/predix_db/GTEx-V6p-HapMap-2016-09-08/TW_%s_0.5.db' % tissue)
            cmd += '--covariance %s ' % ('/data2/users/hclee/data/predix_db/GTEx-V6p-HapMap-2016-09-08/TW_%s.txt.gz' % tissue)
            cmd += '--gwas_folder %s ' % (fn0)
            cmd += '--gwas_file_pattern ".*gz" '
            cmd += '--snp_column SNP --effect_allele_column A1 ' 
            cmd += '--non_effect_allele_column A2 ' 
            cmd += '--beta_column BETA ' 
            cmd += '--pvalue_column P ' 
            cmd += '--se_column SE ' 
            cmd += '--output_file %s' % (fn1)
            cmds.append(cmd)
        
    p.map(run_metaxcan, cmds)
