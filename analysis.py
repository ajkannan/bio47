# Ideas: levenschtein match kmer features with partial points
# Remove from the dataset sequences that are too similar

import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale, normalize
from sklearn.cross_validation import train_test_split

from sklearn import datasets
iris = datasets.load_iris()
irisX = iris.data  # we only take the first two features.
irisY = iris.target

def run_pca(X, dims = 2):
	labels = np.array([str(i) for i in xrange(X.shape[0])])
	Xtrain, Xtest, train_labels, test_labels = train_test_split(X, labels, test_size = 1000, train_size = 3000)
	Xnorm = normalize(Xtrain, axis = 1)
	Xnorm = scale(Xnorm, axis = 0)
	pca = PCA(n_components = dims)
	X_proj = pca.fit_transform(Xnorm)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.scatter(X_proj[:,0], X_proj[:,1])
	for i in xrange(X_proj.shape[0]):
		ax.annotate(train_labels[i], (X_proj[i,0], X_proj[i,1]))
	plt.xlabel("PCA Dimension 1")
	plt.ylabel("PCA Dimension 2")
	plt.title("Cluster Analysis")
	plt.show()
	print "Explained_variance ratio: ", pca.explained_variance_ratio_

def kmer_features(data, k = 9):
	# Threshold kmers based on frequency
	freq_thresh = 7000
	kmer_freqs = {}
	for rna in data:
		for i in xrange(len(rna) - k + 1):
			kmer = rna[i:i + k]
			if kmer not in kmer_freqs:
				kmer_freqs[kmer] = 0.0
			kmer_freqs[kmer] += 1.0

	order = {} # key, position in feature vector
	for key in kmer_freqs:
		if kmer_freqs[key] > freq_thresh:
			order[key] = len(order)
	print len(order)
	X = np.zeros((len(data), len(order)))
	for i in xrange(len(data)):
		rna = data[i]
		for i in xrange(len(rna) - k + 1):
			kmer = rna[i:i + k]
			if kmer in order:
				X[i, order[kmer]] += 1.0
	return X

def kmer_barchart(data, k = 9, remove_low_complexity = True):
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
	for i, rna in enumerate(data):
		if 'N' in rna:
			continue # Occurs 39 times in high confidence dataset
		curr_seq = ''
		start_i = 0
		for lcr in indices[i]:
			curr_seq += rna[start_i:lcr[0]]
			start_i = lcr[1] + 1
			
		curr_seq += rna[start_i:]
		pruned_seqs.append(curr_seq)

	return pruned_seqs

if __name__ == "__main__":
	#kmer_barchart(read_data())
	X = kmer_features(read_data())
	run_pca(X)
	#run_pca(irisX)
