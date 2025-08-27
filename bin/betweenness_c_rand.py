from __future__ import division
import os
import networkx
from networkx.algorithms.centrality.betweenness_subset import \
     _rescale, _accumulate_subset
from networkx.algorithms.centrality.betweenness import \
    _single_source_dijkstra_path_basic as dijkstra_b
from networkx.algorithms.centrality.betweenness import \
    _single_source_shortest_path_basic as shortest_path_b

import copy
import random
import json
import argparse


def create_graph(infile, node_weights=None, edge_cutoff=ec):
    g = networkx.Graph()
    if(infile.endswith(".dat")):
        with open(infile, "rU") as datfile:
            for line in datfile:
                n1, n2, w = line.strip().split("\t")
                w = float(w)

		if w >= edge_cutoff:
                	g.add_edge(n1, n2, weight=1/float(w))

    print('graph built')
    return g


def betweenness_centrality_subsets(G, src_trg,
                                  normalized=False,
                                  weight=None):
    """Compute betweenness centrality for (sources, targets) tuples of nodes.
    Adapted from betweenness_centrality_subset.
    """
    for src, trg in src_trg:
        b = dict.fromkeys(G, 0.0)    # b[v]=0 for v in G
        for s in src:
            # single source shortest paths
            if weight is None:  # use BFS
                S, P, sigma = shortest_path_b(G, s)
            else:  # use Dijkstra's algorithm
                S, P, sigma = dijkstra_b(G, s, weight)
            b = _accumulate_subset(b, S, P, sigma, s, trg)
    return b


def main(args):
	geneset1 = args.geneset1
	geneset2 = args.geneset2
	outname = args.outname
	graph = args.graph
	graph_bg = args.graph_bg
	cur_dir = args.cur_dir
	seedi = args.seedi
	ec = 0
	
	G=create_graph(graph, node_weights=None, edge_cutoff=ec)
	with open(graph_bg) as f_bg:
                list_bg = f_bg.read().splitlines()

	list1=[]
	f1=open(geneset2)
        for line in f1:
                line=line.strip()
		if line in set(list_bg):
			list1.append(line)
	list0=[]
        for s in list1:
		S, P, sigma= dijkstra_b(G, s, weight='weight')
		list0.append([S,P,sigma])

	listv=[]
	f=open(cur_dir+ outname+ '.txt')
	for line in f:
		if line.startswith('#') is False:
			entry=line.strip().split('\t')
			listv.append(entry[0])

	dicc={}
        for key in listv:
                dicc[key]=[]
	fgcnv = open(cur_dir+geneset1[:-4] + '_rand'+ seedi)
	ni = 0
	for line in fgcnv:
		if line != '\n':
			ni +=1
			if ni <= 2001:
				#listt=copy.deepcopy(list0) #deep copy list0 to listt
				listt=[];
				for i in range(0,len(list0)):
					S=copy.deepcopy(list0[i][0])
					listt.append([S,list0[i][1],list0[i][2]])

				list2n=[]
				b = dict.fromkeys(G, 0.0)

				entry = line.strip().split('\t')
				for i in range(0,len(entry)):
					list2n.append(entry[i])
				for i in range(0,len(list0)):
					b = _accumulate_subset(b, listt[i][0],listt[i][1], listt[i][2], list1[i], list2n)
				for key in listv:
					dicc[key].append(b[key])

	
	with open(cur_dir+'intermediate/'+ outname+'_rand%s.json'%seedi,'wb') as fp:
		json.dump(dicc,fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' escription');
    parser.add_argument('--geneset1', type=str, help='geneset1');
    parser.add_argument('--geneset2', type=str, help='geneset2');
    parser.add_argument('--outname', type=str, help='outname', default='result');
    parser.add_argument('--graph', type=str, help='graph');
    parser.add_argument('--graph_bg', type=str, help='graph_bg', default='');
    parser.add_argument('--cur_dir', type=str, help='cur_dir', default='./');
    parser.add_argument('--seedi', type=str, help='random seed', default='1');

    args = parser.parse_args();
    main(args);

