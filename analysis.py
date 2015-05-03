# Ideas: levenschtein match kmer features with partial points
# Remove from the dataset sequences that are too similar

import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
import subprocess
import sys
import random
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale, normalize
from sklearn.cross_validation import train_test_split

def run_pca(X, dims = 2):
	labels = np.array([str(i) for i in xrange(X.shape[0])])
	train_labels = [i for i in xrange(X.shape[0])]
	Xtrain = X
	if X.shape[0] > 4000:
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

def kmer_features(data, klow = 4, khigh = 7):
	# Threshold kmers based on frequency
	kmer_freqs = {}
	for rna in data:
		for k in xrange(klow, khigh + 1):
			for i in xrange(len(rna) - k + 1):
				kmer = rna[i:i + k]
				if kmer not in kmer_freqs:
					kmer_freqs[kmer] = 0.0
				kmer_freqs[kmer] += 1.0

	order = {} # key, position in feature vector
	freq_thresh = 0
	for key in kmer_freqs:
		if kmer_freqs[key] > freq_thresh:
			order[key] = len(order)
	print 'Num kmer features: ', len(order)
	
	X = np.zeros((len(data), len(order)))
	for j in xrange(len(data)):
		rna = data[j]
		for i in xrange(len(rna) - k + 1):
			kmer = rna[i:i + k]
			if kmer in order:
				X[j, order[kmer]] += 1.0
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
	plt.savefig("kmer_barchart.pdf")

def remove_similar_seqs(limit = 500):
	data, names = read_data('./Data/TEMP_lncipedia_3_1_hc.fasta')

	if limit is not None:
		rand_seqs = set()
		while len(rand_seqs) < limit:
			curr_rand = random.randint(0,len(names))
			if curr_rand not in rand_seqs:
				rand_seqs.add(curr_rand)
		temp_data = []
		temp_names = []
		for i in rand_seqs:
			temp_data.append(data[i])
			temp_names.append(names[i])
		data = temp_data
		names = temp_names

	pruned_seqs = {}
	removed_seqs = set()
	similarity_threshold = 50 # Sequences over this fraction identity are not inluded in the dataset
	count = 0
	for name, seq in zip(names, data):
		if name in removed_seqs:
			continue
		# Blast sequence, and then check identity.
		with open('./Data/temp_query.fasta', 'w') as f:
			f.write(name + '\n' + seq + '\n')
		output = subprocess.check_output("blastn -evalue 1000 -db ./Data/lncipediaDb -query ./Data/temp_query.fasta".split())
		sim_seqs_raw = output.split('>')[1:]
		#sims = set()
		for sim_seq in sim_seqs_raw:
			sim_name = '>' + sim_seq[:sim_seq.index('\n')]
			id_start_i = sim_seq.index('Identities')
			id_end_i = sim_seq.index('Gaps')
			identity = int(sim_seq[id_start_i:id_end_i].split('(')[1].split('%')[0])
			if identity > similarity_threshold:
				#sims.add(sim_name)
				removed_seqs.add(sim_name)
		#include = True
		#for sim_name in sims:
		#	if sim_name in pruned_seqs:
		#		include = False
		#if include:
		pruned_seqs[name] = seq
		count += 1
		if count % 100 == 0:
			print 'Finished ' + str(count) + ' BLASTS'

	with open('./Data/CLEANED_lncipedia_3_1_hc.fasta', 'w') as f:
		for seq_name in pruned_seqs:
			f.write(seq_name + '\n' + pruned_seqs[seq_name] + '\n')

	print len(pruned_seqs)
	
def filter_low_complexity():
	data, names = read_data('./Data/lncipedia_3_1_hc.fasta')
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

	# Write pruned data to file to use as database for blast search
	with open('./Data/TEMP_lncipedia_3_1_hc.fasta', 'w') as f:
		for name, seq in zip(names, pruned_seqs):
			f.write(name + '\n' + seq + '\n')

def read_data(filename, limit = None):
	data_raw = open(filename, 'r').readlines()
	data = []
	curr_line = ''
	names = []
	for line in data_raw:
		if len(names) == 0:
			names.append(line.strip())
			curr_line = ''
		elif '>' in line and len(curr_line) != 0:
			data.append(curr_line)
			curr_line = ''
			names.append(line.strip())
		else:
			curr_line += line.strip()
	data.append(curr_line)

	return data, names

if __name__ == "__main__":
	job = int(sys.argv[1])

	if job == 1:
		# Write cleaned database
		filter_low_complexity()
		remove_similar_seqs()

	elif job == 2:
		# Cluster
		seqs, names = read_data('./Data/CLEANED_lncipedia_3_1_hc.fasta')
		X = kmer_features(seqs)
		run_pca(X)

	elif job == 3:
		# Clustering sanity test, check with iris
		from sklearn import datasets
		iris = datasets.load_iris()
		irisX = iris.data  # we only take the first two features.
		irisY = iris.target
		run_pca(irisX)
