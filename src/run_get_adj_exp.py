import os
from multiprocessing import Pool
path = '../processed_data/transposed_expression/'
fn_exp = '../processed_data/transposed_expression/%s-expression.txt' 
fn_covariate = '../processed_data/transposed_covariate/%s-covariate.txt'
fn_save = '../processed_data/adj_expression/%s-adj_expression.txt'

id2tissue = [f.split('-')[0] for f in os.listdir(path) if f[-4:] == '.txt']
#id2tissue = ['test']
pool = Pool(5)

if __name__ == "__main__":
    cmds = []
    for t in id2tissue:
        f0, f1, f2 = fn_exp % t, fn_covariate % t, fn_save % t

        if os.path.isfile(f2):
            continue

        cmd = ' '.join(['Rscript get_adj_exp.R', f0, f1, f2, '\n'])
        #cmds.append(cmd)
        print cmd
    #pool.map(os.system, cmds)