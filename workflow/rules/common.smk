import pandas as pd
import numpy as np


##### load config and sample sheets #####

configfile: "config/config.yaml"

samples = [i.strip() for i in open(config["sample"], 'r')]
chrom_list = [i.strip() for i in open(config["chrom_list"], 'r')]
parity = [i.strip() for i in open(config["parity"], 'r')]
parity_0_or_1 = [i for i in parity if i != "parityall"]
allgwasoutputs = [i.strip() for i in open(config["allgwasoutputs"], 'r')]
motherorchild = [i.strip() for i in open(config["motherorchild"], 'r')]
h2samples = [i for i in allgwasoutputs if i != "effectparity_parityall" and  i != "effectSNPwithInteraction_parityall"]
GE = [i.strip() for i in open(config["features"], 'r')]
pheno_all = [i.strip() for i in open(config["phenotypes"], 'r')]
pheno_main = [i for i in pheno_all if i != "fertility"]
pheno_sup = [i for i in pheno_all if i != "gd" and  i != "bw_zscore"]

def format_duos(d, x):
        print(d.head(10))
        print(d.shape)
        d= d.loc[d["IID"].isin(x.IID.values), :]
        print(d.shape)
        d.drop_duplicates(subset= ["IID"], inplace= True, keep= 'first')
        print(d.shape)
        d= d[["IID"]]
        print(d.shape)
        return d

