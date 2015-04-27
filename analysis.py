import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale, normalize

# Tests
from sklearn import datasets
iris = datasets.load_iris()
Xtest = iris.data


def run_pca(X, dims = 2):
	Xnorm = scale(X, axis = 0)
	pca = PCA(n_components = dims)
	X_proj = pca.fit_transform(Xnorm)
	plt.scatter(X_proj[:,0], X_proj[:,1])
	plt.show()
	print "Explained_variance ratio: ", pca.explained_variance_ratio_

def kmer_features(data, k = 9):
	order_and_freq = {} # key, (position in feature vector, frequency)
	for rna in data:
		for i in xrange(len(rna) - k + 1):
			kmer = rna[i:i + k]
			if kmer not in ordering:
				order_and_freq[kmer] = [len(order_and_freq), 1]
			else:
				order_and_freq[kmer][1] += 1.0

	print ordering[:15]

def kmer_barchart(data, k = 15, remove_low_complexity = True):
	freqs = {}
	for rna in data:
		for i in xrange(len(rna) - k + 1):
			kmer = rna[i:i + k]
			if kmer not in freqs:
				freqs[kmer] = 0.0
			freqs[kmer] += 1.0

	condensed_list = []
	for key in freqs:
		if freqs[key] > 3000:
			condensed_list.append((key, freqs[key]))

	print condensed_list
	fig = plt.figure()
	width = .35
	ind = np.arange(len(condensed_list))

	plt.bar(ind, [kmer[1] for kmer in condensed_list])
	plt.xticks(ind + width / 2, [kmer[0] for kmer in condensed_list])
	
	fig.autofmt_xdate()
	#plt.show()
	plt.savefig("kmer_barchart.pdf")

def read_data():
	data_raw = open('./Data/lncipedia_3_1_hc.fasta', 'r').readlines()
	data = []
	curr_line = ''
	for line in data_raw:
		if '>' in line and len(curr_line) != 0:
			data.append(curr_line)
			curr_line = ''
		else:
			curr_line += line.strip()
			
	lcr_raw = open('./Data/lcr.txt', 'r').readlines() # Low complexity regions
	indices = []
	for line in lcr_raw:
		if '>' in line:
			indices.append([])
		else:
			indices[-1].append([int(i.strip()) for i in line.split('-')])

	pruned_seqs = []
	count = 0
	for i, rna in enumerate(data):
		if 'N' in rna:
			count += 1
			continue
		curr_seq = ''
		start_i = 0
		for lcr in indices[i]:
			curr_seq += rna[start_i:lcr[0]]
			start_i = lcr[1] + 1
			
		curr_seq += rna[start_i:]
		pruned_seqs.append(curr_seq)

	print count
	return pruned_seqs

if __name__ == "__main__":
	kmer_barchart(read_data())
	#run_pca(Xtest)
	#kmer_features(read_data())
