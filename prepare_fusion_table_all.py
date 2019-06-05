import os
import re
import pandas as pd

f=open("list", 'r')
line=[line.rstrip('\n') for line in f]

path=[path+'/final-list_candidate-fusion-genes.txt' for path in line]
fname=[os.getcwd()+'/'+fname for fname in path]	

df_f1 = pd.DataFrame()
df_f2 = pd.DataFrame()

for fn in fname:
	df = pd.read_csv(fn, sep='\t')
	df = df.drop_duplicates(subset=df.columns[0:2])

	df_rec = df[df['Fusion_description'].str.contains('reciprocal', regex=False, na=False)].reset_index(drop=True)
	df_uni = df[~df['Fusion_description'].str.contains('reciprocal', regex=False, na=False)].reset_index(drop=True)

	df_rec = df_rec.iloc[:,0:2]
	df_uni = df_uni.iloc[:,0:2]

	df_f1 = pd.concat([df_f1,df_rec], axis=1)
	df_f2 = pd.concat([df_f2,df_uni], axis=1)

df_f1.to_csv("all_fusion_genes_reciprocal.csv")
df_f2.to_csv("all_fusion_genes_non_reci.csv")
