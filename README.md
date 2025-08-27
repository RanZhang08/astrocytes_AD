# astrocytes_AD


## 1. installation of open source packages

This code was tested on a Linux system in Python 2.7.3. 
Two open source packages need to be installed:
Sleipnir: https://github.com/FunctionLab/sleipnir 
NetworkX: https://networkx.org/ (version 1.11)


## 2. download the astrocyte network
```
cd data/
wget https://astrocyte.princeton.edu/media/networks/astrocyte.dab.gz
gunzip astrocyte.dab.gz
```


## 3. run NetDIFF
### 3.1. generate gold-standards
```
outname=astro_H_AD
outdir=./netdiff/
mkdir -p $outdir
mkdir -p ${outdir}/gstds/
mkdir -p ${outdir}/results/
```

contamination_input: this file contains a list of genes that are known contaminations (e.g. gene markers from unrelated cell types, ribosomal genes in TRAP experiments), will be removed from astrocyte marker genes.
```
contamination_input=./data/contamination_genes.txt
```

diff_input: differential gene expression between case and control, it contains two columns -- gene (human entrez id) and differential expression p-value. If the experiment is not done in human, we need to map genes from mouse to human through orghologs mapping first.
```
diff_input=./data/AD_diff_mat.txt
```

Rscript ./bin/netdiff_generate_gstd.R $outdir $outname $diff_input $contamination_input generate_gstd

### 3.2. run SVM and summarize results
Run SVM using SVMperfer from Sleipnir with five-fold cross-validation.
Features are each gene's network connectivity with all genes in the genome, gold-standards are generated in the previous step. Example gold-standards file: ./data/astro_H_AD_gstd_example.txt.
```
network=./data/astrocyte.dab

for seedi in {1..50}
do
    SVMperfer -l ${outdir}/gstds/${outname}_${seedi}.txt -o ${outdir}/results/${outname}_${seedi}.txt -i $network -a -c 5 -e 10 -t 0.1
done
```

summarize results and generate final NetDIFF prediction
```
Rscript ./bin/generate_NetDiff_gstd.R $outdir $outname $diff_input $contamination_input summarize
```

NetDIFF expected output: a txt file with genome-wide genes and a NetDIFF score predicted. NetDIFF expected run time on a desktop computer: several hours


## 4. run NetPATH
### 4.1. generate betweenness centrality between regional astrocyte signatures and AD genes
example genesets
geneset1: astrocyte-expressed genes in hippocampus
```
geneset1=./data/astro_H_exp.txt
```

geneset2: OMIM and HGMD annotated Alzheimer's disease genes that are expresed in astrocytes
```
geneset2=./data/AD_exp.txt
```

geneset1_bg: the background set of genes that are expressed in astrocytes (mapped from mouse to human) and overlaps with the gene nodes of the astrocyte network
```
geneset1_bg=./data/astro_H_exp_bg.txt
```

geneset2_bg: the background set of genes that are annotated in OMIM and HGMD database
```
geneset2_bg=./data/AD_exp_bg.txt
```

assign output name and diretory
```
outname=astro_H_AD
outdir=./netpath/
mkdir -p $outdir
mkdir -p ${outdir}/intermediate/
mkdir -p ${outdir}/results/
```

convert network from dab to dat file format, and only filter for edges with weights within around top 1 percentile (~0.03 edge weight), using Dat2Dab function in Sleipnir
```
Dat2Dab -i $network -c 0.03 -o ./data/astrocyte_top.dat
```

get background genes in the network
```
Dat2Dab -i ./data/astrocyte_top.dat -E > ./data/astrocyte_top_genes.txt
```

calculate betweenness centrality for genes on shortest path between geneset1 and geneset2, and output the exact shortest paths
```
python ./bin/betweenness_c.py --geneset1 $geneset1 --geneset2 $geneset2 --outname $outname --graph ./data/astrocyte_top.dat --graph_bg ./data/astrocyte_top_genes.txt --cur_dir $outdir
```

### 4.2. run permutation test and get significant intermediate genes
generate permutation geneset for 10000 times, geneset1 is sampled from geneset1_bg (consisted of all genes with similar properties with geneset1), and geneset2 is sampled from geneset2_bg
```
Rscript ./bin/netpath_generate_permutation_geneset.R $geneset2 $geneset2_bg
```

calculate betweenness centrality for permutation genesets, each seed contains 2000 permutations
```
for seedi in 0 1 2 3 4
do
	while [ ! -f ${outdir}/${outname}.txt ]
	do
		sleep 10m
	done

	if [ ! -f ${outdir}/intermediate/${outname}_rand${seedi}.json ]
	then
		python ./bin/betweenness_c_rand.py --geneset1 ${geneset1} --geneset2 ${geneset2//.txt/}_rand${seedi} --outname $outname --graph ./data/astrocyte_top.dat --graph_bg ./data/astrocytes_genes.txt --cur_dir $outdir --seedi $seedi
	fi
done
```

calculate p-value and return genes having statistically significant connection between the two gene sets tested (e.g. hippocampus astrocytes genes and AD genes)
```
python ./bin/netpath_calc_pval.py $outname $outdir
```

NetPATH expected output: a txt file with intermediate genes falling on shortest paths between geneset1 and geneset2, with p-value and FDR calculated (based on permutation test)
NetPATH expected run time on a desktop computer: a day

