from __future__ import division
import os
import sys
import argparse
import networkx

from networkx.algorithms.centrality.betweenness_subset import \
     _rescale, _accumulate_subset
from networkx.algorithms.centrality.betweenness import \
    _single_source_dijkstra_path_basic as dijkstra_b
from networkx.algorithms.centrality.betweenness import \
    _single_source_shortest_path_basic as shortest_path_b


# readdat.py
def create_graph(infile, node_weights=None, edge_cutoff=ec):
	g = networkx.Graph()
	if(infile.endswith(".dat")):
		with open(infile, "rU") as datfile:
			for line in datfile:
				n1, n2, w = line.strip().split("\t")
				w = float(w)
				if node_weights != None:
					if n1 in node_weights and n2 in node_weights:
						adj_weight = w * (node_weights[n1] + node_weights[n2]) / 2
						if adj_weight > 0:
							g.add_edge(n1, n2, weight=adj_weight)
				elif w >= edge_cutoff:
					g.add_edge(n1, n2, weight=1/float(w)) # here take node weight as 1/functional_linkage
	else:
		network = dat.dat(infile)
		for gene in network.gene_list :
			for neighbour in network.get_neighbors(gene, edge_cutoff):
				if not g.has_edge(gene, neighbour):
					g.add_edge(gene, neighbour, weight=network.get_value(network.get_index(gene), network.get_index(neighbour)))

	print('graph built')
	return g


def betweenness_centrality_subsets(G, src_trg,
                                  normalized=False,
                                  weight=None):
	"""Compute betweenness centrality for (sources, targets) tuples of nodes.
	Adapted from betweenness_centrality_subset.
	"""
	# validate
	for src, trg in src_trg:
		# memoize shortest paths computation
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
	ec = 0
	if os.path.isfile(cur_dir + outname +'.txt') is False:
		G=create_graph("%s"%graph, node_weights=None, edge_cutoff=ec)
		with open(graph_bg) as f_bg:
			list_bg = f_bg.read().splitlines()

		with open(geneset1) as f1:
			list1_ori = f1.read().splitlines()

		with open(geneset2) as f2:
			list2_ori = f2.read().splitlines()

		list1 = list(set(list1_ori) & set(list_bg))
		list2 = list(set(list2_ori) & set(list_bg))

		src_trg = [(list1, list2)]
		dict=betweenness_centrality_subsets(G, src_trg, weight="weight")
		
		## output 
		fout = open(cur_dir +outname +'.txt','w')
		fout.write('#gene\tscore\n')
		for key in dict.keys():
			if float(dict[key])!=0:
				fout.write(key+'\t'+str(dict[key])+'\n')
		fout.close()

		## output individual shortest paths
		if True:
			fouti = open(cur_dir + outname+'_indiv.txt','w')
			for g1 in list1:
				for g2 in list2:
					fouti.write(g1+'\t'+g2)
					src_trg=[([g1],[g2])]
					dict=betweenness_centrality_subsets(G, src_trg, weight="weight")
					for key in dict.keys():
						if float(dict[key])!=0:
						       fouti.write('\t'+key)
					fouti.write('\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' escription');
    parser.add_argument('--geneset1', type=str, help='geneset1');
    parser.add_argument('--geneset2', type=str, help='geneset2');
    parser.add_argument('--outname', type=str, help='outname', default='result');
    parser.add_argument('--graph', type=str, help='graph');
    parser.add_argument('--graph_bg', type=str, help='graph_bg');
    parser.add_argument('--cur_dir', type=str, help='cur_dir', default='./');

    args = parser.parse_args();
    main(args);

