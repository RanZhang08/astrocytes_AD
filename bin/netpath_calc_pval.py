import os
import os.path
from scipy import stats
import json
import numpy as np
from sys import argv

outname = argv[1]
indir = argv[2]

cur_dir = indir+'/intermediate/'
out_dir = indir+'/results/'

## read in genes falling within shortest path between the original genesets
dict={}
f=open(out_dir+ outname+ '.txt')
for line in f:
	if line.startswith('#') is False:
		entry=line.split('\t')
		dict[entry[0]]=[]
		dict[entry[0]].append(float(entry[1]))


## compare with betweenness centrality values between randomized genesets and calculate p-value and FDR
if os.path.isfile(cur_dir+ outname+ '_rand1.json') is True: 
	for i in range(0,5):
		if os.path.isfile(cur_dir+ outname+ '_rand'+ str(i)+ '.json') is True: 
			with open(cur_dir+ outname+ '_rand%s.json'%str(i),'rb') as fp:
				dic=json.load(fp)
			
			for key in dic.keys():
				for item in dic[key]:
					dict[key].append(float(item))

	fout=open(cur_dir+ outname+ '_rand10000_pval.txt', 'w')

	for key in dict.keys():
		#N=10000
		N=len(dict[key])
		n=0
		for i in range(0,N):
			if dict[key][i]>=dict[key][0]:
				n+=1
		npval=float(n)/N
		fout.write(str(key)+'\t'+str(dict[key][0])+'\t'+str(npval)+'\n')

	fout.close()

	
	os.system('sort -g -k3 '+cur_dir+ outname+ '_rand10000_pval.txt >out1')

	n=os.popen('wc -l out1').read().split(' ')[0]
	f=open('out1')
	fout=open(out_dir+ outname+ '_rand10000_pval_qval.txt', 'w')
	t=1
	for line in f:
		line=line.strip()
		entry=line.split('\t')
		qval=float(entry[2])*int(n)/float(t)
		fout.write(line+'\t'+str(qval)+'\n')
		t+=1

	fout.close()
	os.system('rm -rf out1')



