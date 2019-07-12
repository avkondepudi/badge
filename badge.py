'''
Badge class
'''

import sys
import numpy as np; np.random.seed(2)
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums, combine_pvalues
from utils import random_color

SHOW_AX_YAXIS_AXES = False
FIG_WIDTH = 16
FIG_HEIGHT = 9

class Badge(object):

	def __init__(self, _input):
		'''
		:param _input: input file name (str) or input as pd.Dataframe()
		'''

		if isinstance(_input, str): _df = pd.read_csv(_input, dtype={'gene': str, 'barcodename': str, 'value': np.float64, 'sample_id': str, 'cell_type': str})
		else: _df = _input

		_df = _df.loc[:, ~_df.columns.str.contains('^Unnamed')]

		self._genes = list(sorted(set(_df['gene'].values)))
		self._samples = list(sorted(set(_df['sample_id'].values)))
		self._ctypes = list(sorted(set(_df['cell_type'].values)))
		self._df = _df.copy()
		del _df

	@staticmethod
	def genes_by_sample(df, genes, ctypes, sample_id, key_vals, return_pd=False):
		'''
		:param df: pd.Dataframe() with data
		:param key_vals: singular, cell type to check for differentially expressed genes
		
		:return marray/mdf: candidate genes by sample
		'''

		assert (key_vals==list(set(key_vals).intersection(ctypes))), "%s cell type not found in data." % (key_vals[0])
	
		marray = np.zeros((len(genes), len([c for c in ctypes if c not in key_vals])))
		
		_df = df[(df['sample_id']==sample_id)]
		for i, gene in enumerate(genes):
			gdf = _df[_df['gene']==gene]
			for key_val in key_vals:
				
				pvals = []
				
				for c in ctypes:
					if c not in key_vals:
						vals = ranksums(gdf[gdf['cell_type']==key_val]['value'].dropna().values, gdf[gdf['cell_type']==c]['value'].dropna().values)
						tval, pval = vals[0], vals[1]
						if tval <= 0: pvals.append(1-0.5*pval)
						else: pvals.append(0.5*pval)
				
				for j in range(len(pvals)): marray[i][j]=pvals[j]		
		
		if return_pd:
			
			mdf = pd.DataFrame(columns=[c for c in ctypes if c not in key_vals])
			
			i=0
			for c in ctypes:
				if c not in key_vals:
					mpd[c] = marray[:,i]
					i+=1

			mdf['gene'] = genes
			mdf.set_index('gene', inplace=True)

			return mdf
		
		else: return marray

	@staticmethod
	def genes(df, genes, ctypes, samples, key_vals):
		'''
		:param key_vals: singular, cell type to check for differentially expressed genes

		:return: candidate genes
		'''
	
		marray = np.zeros((len(genes), len([c for c in ctypes if c not in key_vals]), len(samples)))
		
		for i, sample in enumerate(samples):
			marray[:, :, i] = Badge.genes_by_sample(df, genes, ctypes, sample, key_vals)
		
		narray = np.zeros((len(genes), len([c for c in ctypes if c not in key_vals])))
		
		for i, gene in enumerate(genes):
			for j, c in enumerate([u for u in ctypes if u not in key_vals]):
				narray[i][j] = combine_pvalues(marray[i][j][~np.isnan(marray[i][j])])[1]

		gl = []
		q = 0.05
		genes = np.array(genes)
		for i in range(len([u for u in ctypes if u not in key_vals])):
			pvals = narray[:, i]
			sorted_pvals = pvals[pvals.argsort()]
			sorted_gnames = genes[pvals.argsort()]

			j = np.arange(1, len(sorted_pvals)+1)

			below = sorted_pvals < (q * j / len(sorted_pvals))
			try: 
				max_below = np.max(np.where(below)[0])
				gl.append(sorted_gnames[:max_below+1])
			except ValueError:
				gl.append([""])
			
		return list(set(gl[0]).intersection(*gl))

	def candidate_genes(self, ctypes):
		'''
		candidate genes
		'''

		mdict = {}
		for ctype in ctypes:
			mdict[ctype] = self.genes(self._df, self._genes, self._ctypes, self._samples, [ctype])
		return mdict

	def graph(self, genes, save=False):
		'''
		:return None: graphs expression of genes across samples
		'''

		if isinstance(genes, str):
			_genes = []
			_genes.append(genes)
			genes = _genes

		for gene in genes:

			df = self._df[(self._df['gene']==gene)]

			fig = plt.figure(gene, figsize=(FIG_WIDTH, FIG_HEIGHT))
			ax = fig.add_subplot()

			for c in self._ctypes: ax.scatter([df[(df['cell_type']==c) & (df['sample_id']==s)]['value'].mean() for s in self._samples], [range(len(self._samples))[u]+np.random.random_sample()/5. for u in range(len(self._samples))], s=100, color=random_color(), label=c, alpha=0.5)

			yticloc = []
			bars = ax.barh(range(len(self._samples)), np.ones(len(self._samples)), 0)
			[yticloc.append(bar.get_y()) for bar in bars]

			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.tick_params(axis='both', which='both', length=5)
			ax.set_xlim([0,max(df['value'].values)*(1.01)])
			ax.set_yticks(yticloc)
			ax.set_yticklabels(self._samples)
			ax.yaxis.grid(SHOW_AX_YAXIS_AXES)
			ax.set_xlabel('EXPRESSION')
			ax.set_ylabel('SAMPLE')
			ax.legend(title='CELL TYPE')
			ax.legend(title='CELL TYPE', loc='center left', bbox_to_anchor=(1, 0.5))        

			plt.suptitle(gene, fontweight='bold')

			if not save:
				while True:
					try:
						plt.show()
					except UnicodeDecodeError:
						continue
					break
			else: plt.savefig('./figs/'+str(gene))

			del fig, ax

		return None
