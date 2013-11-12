from antz.io import *
import numpy as np
import matplotlib.pyplot as plt
import itertools as it


def isSeq( seq, alphabet ):
	for char in seq.upper():
		if char not in alphabet and char is not 'N':
			return False
	return True 

def seqType( seq ):
	seqs = { 'D':1, 'R':1, 'P':1 }
	for char in seq.upper():
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

def isRNA( seq ):
	return isSeq( seq, alpabet=RNA.alphabet )

def isDNA( seq ):
	return isSeq( seq, alphabet=DNA.alphabet )

def isProtein( seq ):
	return isSeq( seq, alphabet=Protein.alphabet )

def rc( seq, compliment=None ):
	if not compliment:
		compliment = seqType(seq).compliment
	return map( lambda x: _rc[x], seq.upper() )

def composition( sequence, alphabet=None, window=None ):
	n = len( sequence )
	if not window:
		window = n 
	if not alphabet:
		alphabet = seqType( seq ).alphabet

	for i in xrange( n-window+1 ):
		# If first window, must run through sequence
		# If a window, you want a list, otherwise just an integer
		if i == 0: 
			composition = { i:0 for i in alphabet }
			for char in sequence[ :window ]:
				composition[char] += 1.
			for char in alphabet:
				composition[char] /= window
			if window != n:
				for char in alphabet:
					composition[char] = [ composition[char] ]

		# If not first window, the percentages for only at most two nucleotides
		# change, since you're sliding by one. You can just change those and
		# go through very quickly that way. 	
		else:
			last_out = sequence[i-1]
			next_in = sequence[i+window-1]
			for letter in alphabet:
				freq = composition[letter][-1]
				if letter == last_out:
					freq -= 1./window
				if letter == next_in:
					freq += 1./window
				composition[letter].append( freq )

	return composition 

def clean( self, alphabet=None, replacement='N' ):
	if not alphabet:
		alphabet = seqType( seq ).alphabet
	return map( lambda x: x if x in alphabet else 'N', seq.upper() )

class Sequence( object ):
	def __init__( self, sequence, comments=None, **kwargs ):
		self.sequence = sequence
		self.comments = comments
		self.n = len( self.sequence )
		for key, value in kwargs.iteritems():
			print key, value
			if not hasattr( self, key ):
				setattr( self, key, value )

	def clean( self, alphabet=None, replacement='N' ):
		return clean( self.sequence, alphabet, replacement )
	
	def to_fasta( self, filename, comments=None ):
		FastA(sequence=self.sequence, comments=comments or self.comments or '').to_file(filename)

	def composition( self, alphabet=None, window=None ):
		alphabet = alphabet or self.__class__.alphabet
		return composition( self.sequence, alphabet=alphabet, window=window )

	def plot( self, x, window=None, title=None, c='c', alpha=0.8, linewidth=2, **kwargs ):
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
		add = lambda x, y: x + self.__class__._mW[y]
		if not window:
			return reduce( add, self.sequence, 0 )
		return [ reduce( add, self.sequence[i:i+window], 0 ) for i in xrange(self.n-window+1) ]

	def reverse( self ):
		return self.__class__( self.sequence[::-1] )

	def shuffle( self ):
		self.sequence = list( self.sequence )
		np.random.shuffle( self.sequence )
		self.sequence = ''.join( seq ) 

	def kmers( self, k ):
		if type(k) == tuple:
			m, k = k
			mult=True

		alpha = self.__class__.alphabet
		kmers = { ''.join(kmer):0 for kmer in it.product( alpha, repeat=k )} 
		for i in xrange( self.n-k+1 ):
			kmer = self.sequence[ i:i+k ]
			kmers[ kmer ] += 1

		if mult:
			for i in np.arange( k-1, m-1, -1 ):
				print i
				imers = it.product( alpha, repeat=i )
				for imer in imers:
					imer = ''.join(imer)
					kmers[ imer ] = np.sum([ kmers[ imer+char ] for char in alpha ])

		for key, val in kmers.items():
			print key, val
		return kmers

	@classmethod
	def random( cls, n, alphabet=None ):
		if not alphabet:
			alphabet = cls.alphabet
		alpha = list( cls.alphabet )
		return cls( ''.join([ alpha[i] for i in np.random.randint( len(cls.alphabet), size=n) ]) )

	@classmethod
	def from_fasta( cls, filename ):
		return ( cls( seq ) for _, seq in FastA(seq) )


class DNA( Sequence ):
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
		add = lambda x, y: x + ( y in 'GC' )
		m = self.n-window+1
		if not window:
			return reduce( add, self.sequence, 0. ) / self.n
		return [ reduce( add, self.sequence[i:i+window], 0. )/window for i in xrange(m) ]

	def rc( self, sequence=None ):
		return rc( self.sequence, compliment=self.__class__.compliment )

	def transcribe( self ):
		return RNA( map( lambda x: x if x is not 'T' else 'U', self.sequence ) )

	def translate( self ):
		return Protein( map( lambda x: DNA.codons[x], self.sequence ) )

	def ORFs( self ):
		for i in xrange( self.n-2 ):
			if DNA.codons[ self.sequence[ i:i+3 ] ] == '*':
				for j in xrange( 3, i, 3 ):
					if self.sequence[ i-j:i-j+3 ] == 'ATG':
						seq = self.sequence[ i-j+3:i+3 ]
						yield DNA( seq, start=i-j+3, end=i+3, length=j, frame=i%3 )
						break

			if DNA.codons[ DNA.rc( self.sequence[ i:i+3 ] ) ] == '*':
				for j in xrange( 3, i, 3 ):
					if DNA.rc( self.sequence[ i+j:i+j+3 ] ) == 'ATG':
						seq = self.sequence[ i:i+j+3 ]
						yield DNA( seq, start=self.n-i-j-3, end=self.n-i, length=j+3, frame=-(self.n-i)%3)
						break

class RNA( Sequence ):
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
		return Protein( map( lambda x: RNA.codons[x], self.sequence ) )

	def rc( self ):
		return RNA( rc( self.sequence, compliment=self.__class__.compliment ) )

class Protein( Sequence ):
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

