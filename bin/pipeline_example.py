import os
import pandas as pd

#you can use your own pipeline tools to make a simple pipeline to draw many different chromosome areas at the same time. Here is an example:

df=pd.read_table('BxPC3.bed',sep='\t',names=['chrs','start','end','name'])
for i in range(len(df)):
  cmd='python3 Annoroad_Browser.py -fl show.list -o {0} -chr {1} -st {2} -ed {3} -b1 {4} -b2 {4}'.format(df['name'][i],df['chrs'][i],df['start'][i],df['end'][i],1000000)
  os.system(cmd)
  
