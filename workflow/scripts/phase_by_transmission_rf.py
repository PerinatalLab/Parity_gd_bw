import pandas as pd
import numpy as np
import csv


PREG_ID= 'fam'
Sentrix= 'IID'

def format_df(df):
    df[['chr', 'pos', 'ref', 'eff']]= df['index'].str.split(':', expand= True)
    cols = list(df.columns.values)
    cols= cols[-4:] + cols[:-4]
    df= df[cols]
    df.drop(['index'], axis= 1, inplace= True)
    return df


d= pd.read_csv(snakemake.input[2], sep= '\t', header= 0) #duos
d.dropna(axis= 0, inplace= True)
print(d.shape)

with open(snakemake.input[0], 'r') as infile:
    reader= csv.DictReader(infile, delimiter= '\t')
    fets_cols= reader.fieldnames # columns = child ID

with open(snakemake.input[1], 'r') as infile:
    reader= csv.DictReader(infile, delimiter= '\t')
    moms_cols= reader.fieldnames #columns = mother ID

d= d.loc[d.Child.isin(fets_cols), :]
d= d.loc[d.Mother.isin(moms_cols), :]
print(d.shape)

fets_snp_list= list()
h1_df_list= list()
h3_df_list= list()

for fets in pd.read_csv(snakemake.input[0], sep='\t', header=0, chunksize=100):
    fets_snp_list.append(fets.chr.apply(str) + ':' + fets.pos.apply(str) + ':' + fets.ref + ':' + fets.eff)
    fets = fets[d.Child]
    fets = fets.astype(str)
    fets = np.where(fets == '0', '0|0', np.where(fets == '1', '1|1', fets))
    h3 = np.where((fets == '0|0') | (fets == '0|1'), 0, np.where((fets == '1|0') | (fets == '1|1'), 1, np.nan))
    h1 = np.where((fets == '0|0') | (fets == '1|0'), 0, np.where((fets == '0|1') | (fets == '1|1'), 1, np.nan))
    h1_df_list.append(h1)
    h3_df_list.append(h3)

varnames= pd.concat(fets_snp_list).values.tolist()
h1= np.concatenate(h1_df_list)
h3= np.concatenate(h3_df_list)
print(h1.shape)

moms_df_list= list()

for moms in pd.read_csv(snakemake.input[1], sep= '\t', header= 0, chunksize= 100):
    moms= moms[d.Mother]
    moms= np.where(
            moms == '0|0', 0,
            np.where(
                (moms == '0|1') | (moms == '1|0'), 1,
                np.where(moms == '1|1', 2, np.nan
                         )
                )
            )
    moms_df_list.append(moms)
moms= np.concatenate(moms_df_list)

h2= moms - h1

d['fam'] = d['Child'] + ':' + d['Mother']

h1= pd.DataFrame(data= h1, columns= d[PREG_ID], index= varnames).reset_index()
h2= pd.DataFrame(data= h2, columns= d[PREG_ID], index= varnames).reset_index()
h3= pd.DataFrame(data= h3, columns= d[PREG_ID], index= varnames).reset_index()

h1.to_csv(snakemake.output[0], header= True, sep= '\t', index= False)
h2.to_csv(snakemake.output[1], header= True, sep= '\t', index= False)
h3.to_csv(snakemake.output[2], header= True, sep= '\t', index= False)
