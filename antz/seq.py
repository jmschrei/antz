# seq.py
# Contact: Jacob Schreiber
#          jmschr@cs.washington.edu

'''
Provides the core sequence data types for the other libraries. Handles
sequences of DNA, RNA, and protein, and common computational things you may
want to do with them. These data types will form the core of other modules
such as the alignment module or the io module.
'''

import numpy as np
import matplotlib.pyplot as plt
import itertools as it

def seqType( sequence ):
	'''
	Determine the most likely sequence type for this sequence. It does this by
	going through the sequence and looking for a character which does not belong
	in each sequence type until only one sequence type is left. If there are
	several choices left (for example, 'ACTACAACTACTACATCAC'), then assume that
	it is DNA, then RNA if DNA is not a choice.
	'''

	seqs = { 'D':1, 'R':1, 'P':1 }
	for char in sequence:
		if seqs['D'] and char not in DNA.alphabet:
			seqs['D'] = 0
		elif seqs['R'] and char not in RNA.alphabet:
			seqs['R'] = 0
		elif seqs['P'] and char not in Protein.alphabet:
			seqs['P'] = 0

		if sum( seqs.values() ) == 1:
			break 

	if seqs['D']:
		return DNA
	if seqs['R']:
		return RNA
	if seqs['P']:
		return Protein

def isSeq( seq, alphabet ):
	'''
	Returns whether or not the sequence qualifies under the alphabet given.
	'''

	for char in seq.upper():
		if char not in alphabet and char is not 'N':
			return False
	return True 

def isRNA( sequence ):
	'''
	Return whether or not this sequence qualifies as RNA.
	'''

	return isSeq( sequence, alpabet=RNA.alphabet )

def isDNA( sequence ):
	'''
	Return whether or not this sequence qualifies as DNA.
	'''

	return isSeq( sequence, alphabet=DNA.alphabet )

def isProtein( sequence ):
	'''
	Return whether or not this sequence qualifies as protein.
	'''

	return isSeq( sequence, alphabet=Protein.alphabet )

def rc( seq, compliment=None ):
	'''
	Return the reverse compliment of the strand. Either define a compliment
	mapping using a dictionary, or it will try to find the best given the
	composition of the sequence given.
	'''

	if not compliment:
		compliment = seqType(seq).compliment
	return ''.join( map( lambda x: compliment[x], seq.upper() ) )[::-1]

def composition( sequence, alphabet=None, window=None ):
	'''
	Return the composition of the sequence. This can be done over the entire
	file, or as a sliding window over the sequence. If it is over the entire
	sequence, a dictionary of keys-percentages is returned, otherwise a
	dictionary of key-lists is provided, where the lists are of percentages
	using the sliding window
	'''

	n = len( sequence )
	if not window:
		window = n 
	if not alphabet:
		alphabet = seqType( seq ).alphabet

	# If no window, just return the composition across the entire sequence
	if not window:
		composition = { i : 0. for i in alphabet }

		# Go through each character and add a count
		for char in sequence:
			composition[ char ] += 1

		# Now normalize to get percentages
		for char in alphabet:
			composition[ char ] /= window 

	# If a sliding window, returning composition using this sliding window 
	else:
		composition = { i : [0.] for i in alphabet }
		
		# Get starting probabilities for the first window
		for char in sequence[:window]:
			composition[char] += 1. / window
		
		# Go through the remainder of the sequence
		for i in xrange( n-window+1 ):
			last_out = sequence[i-1]
			next_in  = sequence[i+window-1]

			for letter in alphabet:
				# Get the last probability from the window
				prob = composition[letter][-1]
				if letter == last_out:
					prob -= 1./window
				if letter == next_in:
					prob += 1./window
				composition[letter].append( prob )

	return composition 

def clean( sequence, alphabet=None, replacement='N' ):
	'''
	Clean a sequence by turning all characters not belonging in the alphabet
	into Ns. If no alphabet is provided, then attempt to determine if the type
	is DNA, RNA, or Protein automatically based on the composition of
	characters. 
	'''

	if not alphabet:
		alphabet = seqType( sequence ).alphabet
	return map( lambda x: x if x in alphabet else 'N', seq.upper() )

class Sequence( object ):
	'''
	A generic sequence object. This stores the sequence information and whatever
	other metadata is needed on the object. 
	'''

	def __init__( self, sequence, **kwargs ):
		'''
		Take in the sequence and any other data, and save it.
		'''

		self.sequence = sequence
		self.attributes = kwargs.keys()
		for key, value in kwargs.iteritems():
			if not hasattr( self, key ):
				setattr( self, key, value )

	def __str__( self ):
		'''
		A string representation of the sequence is just a string representation
		of the underlying sequence.
		'''

		return self.to_fasta( self.attributes )

	@property
	def n( self ):
		return len( self.sequence )

	def clean( self, alphabet=None, replacement='N' ):
		'''
		Clean the sequence by replacing all characters not in the alphabet with
		'N's.
		'''

		return clean( self.sequence, alphabet, replacement )
	
	def composition( self, alphabet=None, window=None ):
		'''
		Return the composition of this sequence, either as a full sequence or
		along a sliding window.
		'''

		alphabet = alphabet or self.__class__.alphabet
		return composition( self.sequence, alphabet=alphabet, window=window )

	def plot( self, x, window=None, title=None, c='c', alpha=0.8, linewidth=1.5, 
		**kwargs ):
		'''
		Plot one of many characteristics of this sequence. x should be a string
		of the function name you want to plot. 
		'''

		x, xlabel = getattr( self, x )( window=window ), x

		if type(x) == dict:
			m = len( x.keys() )

			if all([ type( value ) == float for value in x.values() ]):
				plt.figure( facecolor='0.9' )
				plt.subplot( 111, axisbg='0.8' )
				plt.title( xlabel )
				plt.bar( np.arange(m)+0.125, x.values(), color='c', alpha=0.75, width=0.75 )
				plt.xticks( np.arange(m)+0.5, x.keys() )

			else:
				fig = plt.figure( figsize=( 6, 6*min(m, 3) ), facecolor='0.9' )
				fig.text( 0.5, .95, xlabel, horizontalalignment='center', verticalalignment='top' )
				for i, (key, value) in enumerate( x.iteritems() ):
					plt.subplot( m, 1, i, axisbg='0.8' )
					plt.title( key )
					plt.plot( xrange( len(value) ), value, color=c, alpha=alpha, 
															linewidth=linewidth, **kwargs )
					plt.grid( True, alpha=0.75 )
					plt.xlim( [0, self.n-window+1] )

		elif type(x) == list:
			fig = plt.figure( facecolor='0.9'  )
			fig.text( 0.5, 0.95, xlabel, horizontalalignment='center', verticalalignment='top' )
			plt.subplot( 111, axisbg='0.8' )
			plt.plot( xrange(len(x)), x, color=c, alpha=alpha, linewidth=linewidth, **kwargs )
			plt.grid( True, alpha=0.75 )
			plt.xlim([ 0, self.n ])

		plt.show()

	def mW( self, window=None ):
		'''
		Return the molecular weight of the sequence.
		'''

		add = lambda x, y: x + self.__class__._mW[y]
		if not window:
			return reduce( add, self.sequence, 0 )
		return [ reduce( add, self.sequence[i:i+window], 0 ) for i in xrange(self.n-window+1) ]

	def reverse( self ):
		'''
		Return the reverse of the sequence.
		'''

		return self.__class__( self.sequence[::-1] )

	def shuffle( self ):
		self.sequence = list( self.sequence )
		np.random.shuffle( self.sequence )
		self.sequence = ''.join( seq ) 

	def kmers( self, k ):
		'''
		Return all kmers of length k. 
		'''

		alpha = self.__class__.alphabet
		kmers = { ''.join(kmer):0 for kmer in it.product( alpha, repeat=k )} 
		for i in xrange( self.n-k+1 ):
			kmer = self.sequence[ i:i+k ]
			kmers[ kmer ] += 1

		return kmers

	def to_fasta( self, attrs=[] ):
		'''
		Return a FastA entry for this sequence. Give a sequence of attributes
		to be returned in the comment line of the FastA file, or else nothing
		will be.
		'''

		return '>' + ' '.join( getattr( self, attr ) for attr in attrs ) + \
			'\n' + self.sequence  + '\n'

	@classmethod
	def random( cls, n, alphabet=None ):
		'''
		Return a random sequence from an alphabet assuming uniform probabilities
		on each of the characters
		'''

		if not alphabet:
			alphabet = cls.alphabet
		return cls( ''.join([ alphabet[i] for i in np.random.randint( 
			len(cls.alphabet), size=n) ]) )


class DNA( Sequence ):
	'''
	This is specifically a DNA sequence, which means that its alphabet is
	clamped to the four canonical nucleotides, their molecular weights
	and compliments are set, and a list of codons is provided.
	'''

	alphabet = set('ACGT')
	_mW = { 'A': 331.2, 'C': 307.2, 'G': 347.2, 'T': 322.2 }
	compliment = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }
	codons = {  'CTT':'L', 'TAG':'*', 'ACA':'T', 'ACG':'T', 'ATC':'I', 'AAC':'N', 'ATA':'I',
				'AGG':'R', 'CCT':'P', 'ACT':'T', 'AGC':'S', 'AAG':'K', 'AGA':'R', 'CAT':'H', 
				'AAT':'N', 'ATT':'I', 'CTG':'L', 'CTA':'L', 'CTC':'L', 'CAC':'H', 'AAA':'K', 
				'CCG':'P', 'AGT':'S', 'CCA':'P', 'CAA':'Q', 'CCC':'P', 'TAT':'Y', 'GGT':'G', 
				'TGT':'C', 'CGA':'R', 'CAG':'Q', 'TCT':'S', 'GAT':'D', 'CGG':'R', 'TTT':'F', 
				'TGC':'C', 'GGG':'G', 'TGA':'*', 'GGA':'G', 'TGG':'W', 'GGC':'G', 'TAC':'Y', 
				'TTC':'F', 'TCG':'S', 'TTA':'L', 'TTG':'L', 'TCC':'S', 'ACC':'T', 'TAA':'*', 
				'GCA':'A', 'GTA':'V', 'GCC':'A', 'GTC':'V', 'GCG':'A', 'GTG':'V', 'GAG':'E', 
				'GTT':'V', 'GCT':'A', 'GAC':'D', 'CGT':'R', 'GAA':'E', 'TCA':'S', 'ATG':'M', 
				'CGC':'R' }

	def GC( self, window=None ):
		'''
		Return GC content. This can be either a single number representing the
		GC content across the entire sequence, or as a sliding window.
		'''

		add = lambda x, y: x + ( y in 'GC' )
		m = self.n-window+1
		if not window:
			return reduce( add, self.sequence, 0. ) / self.n
		return [ reduce( add, self.sequence[i:i+window], 0. )/window for i in xrange(m) ]

	def rc( self ):
		'''
		Return the reverse compliment of this sequence.
		'''

		return rc( self.sequence, compliment=self.__class__.compliment )

	def transcribe( self ):
		'''
		Return the RNA object which would have been produced if this were transcribed.
		'''

		return RNA( map( lambda x: x if x is not 'T' else 'U', self.sequence ) )

	def translate( self ):
		'''
		Go all the way to protein and convert the sequence to protein using the
		codon table.
		'''

		return Protein( map( lambda x: DNA.codons[x], self.sequence ) )

	def ORFs( self ):
		'''
		Do a simple ORF finding algorithm, yielding ORFs one at a time. Start
		ORF finding at the start codon, and end when it reaches an end codon.
		Then go across the reverse compliment and look for ORFs on that as well.
		'''

		for i in xrange( self.n-2 ):
			if DNA.codons[ self.sequence[ i:i+3 ] ] == '*':
				for j in xrange( 3, i, 3 ):
					if self.sequence[ i-j:i-j+3 ] == 'ATG':
						seq = self.sequence[ i-j:i+3 ]
						yield DNA( seq, start=i-j, end=i+3, length=j, frame=i%3 )
						break

			if DNA.codons[ rc( self.sequence[ i:i+3 ] ) ] == '*':
				for j in xrange( 3, self.n-i, 3 ):
					if rc( self.sequence[ i+j:i+j+3 ] ) == 'ATG':
						seq = rc( self.sequence[ i:i+j+3 ] )
						yield DNA( seq, start=self.n-i-j-3, end=self.n-i, length=j+3, frame=-((self.n-i)%3))
						break

class RNA( Sequence ):
	'''
	This is specifically an RNA sequence, which means that its alphabet is
	clamped to the four canonical nucleotides, their molecular weights
	and compliments are set, and a list of codons is provided.
	'''

	alphabet = set('ACGU')
	_mW = { 'A': 347.2, 'C': 323.2, 'G': 363.2, 'U': 324.2 }
	compliment = { 'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G' }
	codons = { 'ACC':'T', 'GUC':'V', 'ACA':'T', 'ACG':'T', 'GUU':'V', 'AAC':'N', 'CCU':'P',
			   'UAU':'Y', 'AGC':'S', 'AUC':'I', 'CAU':'H', 'AAU':'N', 'AGU':'S', 'ACU':'T', 
			   'GUG':'V', 'CAC':'H', 'AAA':'K', 'CCG':'P', 'CCA':'P', 'CAA':'Q', 'CCC':'P', 
			   'GGU':'G', 'UCU':'S', 'GCG':'A', 'UGC':'C', 'CAG':'Q', 'UGA':'*', 'UGG':'W', 
			   'CGG':'R', 'UCG':'S', 'AGG':'R', 'GGG':'G', 'UCC':'S', 'UCA':'S', 'GAA':'E', 
			   'UAA':'*', 'GGA':'G', 'UAC':'Y', 'CGU':'R', 'UAG':'*', 'AUA':'I', 'GCA':'A', 
			   'CUU':'L', 'GGC':'G', 'AUG':'M', 'CUG':'L', 'GAG':'E', 'CUC':'L', 'AGA':'R', 
			   'CUA':'L', 'GCC':'A', 'AAG':'K', 'GAU':'D', 'UUU':'F', 'GAC':'D', 'GUA':'V', 
			   'CGA':'R', 'GCU':'A', 'UGU':'C', 'AUU':'I', 'UUG':'L', 'UUA':'L', 'CGC':'R', 
			   'UUC':'F' }

	def translate( self ):
		'''
		Translate this strand directly into protein, using the RNA codon table.
		'''

		return Protein( map( lambda x: RNA.codons[x], self.sequence ) )

	def rc( self ):
		'''
		Return the reverse compliment of this strand. 
		'''

		return RNA( rc( self.sequence, compliment=self.__class__.compliment ) )

class Protein( Sequence ):
	'''
	This is specifically a protein sequence, which means that its alphabet is
	clamped to the twenty canonical amino acids, their molecular weights, pkas,
	and hydropathys are set.
	'''

	alphabet = set('ACDEFGHIKLMNPQRSTVWY')

	_mW = {'A': 89.1, 'C': 121.2, 'D': 133.1, 'E': 147.1, 'F': 165.2, 'G': 75.1,
		'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2, 'M': 149.2, 'N': 132.1, 
		'P': 115.1, 'Q': 146.1, 'R': 174.2, 'S': 105.1, 'T': 119.1, 'V': 117.1, 
		'W': 204.2, 'Y': 181.2  }

	_pKa_c = {'G': 2.34, 'A': 2.34, 'V': 2.32, 'L': 2.36, 'I': 2.36, 'M': 2.28,
			 'P': 1.99, 'F': 1.83, 'W': 2.82, 'N': 2.02, 'Q': 2.17, 'S': 2.21,
			 'T': 2.09, 'Y': 2.20, 'C': 1.96, 'D': 1.88, 'E': 2.19, 'L': 2.18,
			 'R': 2.17, 'H': 1.82 }

	_pKa_n = {'G': 9.60, 'A': 9.69, 'V': 9.62, 'L': 9.60, 'I': 9.60, 'M': 9.21,
			 'P': 10.6, 'F': 9.13, 'W': 9.39, 'N': 8.80, 'Q': 9.13, 'S': 9.15,
			 'T': 9.10, 'Y': 9.11, 'C': 8.18, 'D': 9.60, 'E': 9.67, 'L': 8.95,
			 'R': 9.04, 'H': 9.17 }

	_pKa = { 'D': 3.9, 'E': 4.3, 'R': 12.0, 'K': 10.5, 'H': 6.08, 'C': 8.28, 
			'Y': 10.1 }

	# Hydropathy scores found at http://gcat.davidson.edu/DGPB/kd/aminoacidscores.htm
	_hydropathy = { 'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9,
					'A': 1.8, 'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8,
					'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'Q': -3.5,
					'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5 }

	def hydropathy( self, window=None ):
		'''
		Return the hydropathy of either the entire protein, or a set window.
		'''

		add = lambda x, y: x + self.__class__._hydropathy[y]
		if not window:
			window = self.n
		return [ reduce( add, self.sequence[i:i+window], 0 ) for i in xrange( self.n-window+1 ) ]