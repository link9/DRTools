import pandas as pd
import pubchempy as pcp


def pubchemname_to_cid(name):
    cid = 1
    return cid


#%%    
df = pd.read_csv('../CodeMapping/CMAP-Drugs-1Target-TextMining-filtering.txt', sep='\t', encoding='utf-8')

df = df.loc[:,['CID','DRUG']].drop_duplicates()




#%%

dsigdb = pd.read_csv('D:\data.biodb\DSigDB\DSigDB_All_detailed.txt', sep='\t')
pubchemname = dsigdb.loc[dsigdb.Source == 'D1 PubChem','Drug']
pubchemname = list(pubchemname.unique())
pcp.get_cids('Lenvatinib', 'name')
