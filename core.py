import numpy as np

class Sequence( str ):
	def __init__( self, sequence, comments=None, **kwargs ):
		self.sequence = sequence
		self.comments = comments
		self.n = len( self.sequence )
		for key, value in kwargs.iteritems():
			if not hasattr( self, key ):
				setattr( self, key, value )

	def clean( self, alphabet=None, replacement='N' ):
		if not alphabet:
			alphabet = self.__class__.alphabet
		self.sequence = map( lambda x: x if x in alphabet else replacement, self.sequence )
	
	def to_fasta( self, filename, comments=None ):
		try:
			with open( filename, 'w' ) as outfile:
				outfile.write( comments or self.comments + "\n" )
				outfile.write( self.sequence + "\n" )
		except:
			raise IOError( "Cannot open file {}".format( filename ) )

	@property
	def mW( self ):
		return reduce( lambda x, y: x += self.mW[y], sel.sequence, 0 )

	def 

	@classmethod
	def random( self, n, alphabet=None ):
		if not alphabet:
			alphabet = self.alphabet
		alpha = list( self.alphabet )
		return self( ''.join([ alpha[i] for i in np.random.randint( len(self.alphabet), size=n) ]) )

	@classmethod
	def from_fasta( self, filename ):
		return ( self( seq ) for _, seq in FastA.seq )


class DNA( Sequence ):
	alphabet = set('ACGT')
	mW = { 'A': 331.2, 'C': 307.2, 'G': 347.2, 'T': 322.2 }

	def GC( self ):
		return reduce( lambda x, y: x + ( y in 'GC' ), self.sequence, 0. ) / self.n

class RNA( Sequence ):
	alphabet = set('ACGU')	
	mW = { 'A': 347.2, 'C': 323.2, 'G': 363.2, 'U': 324.2 }

class Protein( Sequence ):
	alphabet = set('ACDEFGHIKLMNPQRSTVWY')

	mW = {'A': 89.1, 'C': 121.2, 'D': 133.1, 'E': 147.1, 'F': 165.2, 'G': 75.1,
		'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2, 'M': 149.2, 'N': 132.1, 
		 'P': 115.1, 'Q': 146.1, 'R': 174.2 'S': 105.1, 'T': 119.1, 'V': 117.1, 
		'W': 204.2, 'Y': 181.2  }

	pKa_c = {'G': 2.34, 'A': 2.34, 'V': 2.32, 'L': 2.36, 'I': 2.36, 'M': 2.28,
			 'P': 1.99, 'F': 1.83, 'W': 2.82, 'N': 2.02, 'Q': 2.17, 'S': 2.21,
			 'T': 2.09, 'Y': 2.20, 'C': 1.96, 'D': 1.88, 'E': 2.19, 'L': 2.18,
			 'R': 2.17, 'H': 1.82 }

	pKa_n = {'G': 9.60, 'A': 9.69, 'V': 9.62, 'L': 9.60, 'I': 9.60, 'M': 9.21,
			 'P': 10.6, 'F': 9.13, 'W': 9.39, 'N': 8.80, 'Q': 9.13, 'S': 9.15,
			 'T': 9.10, 'Y': 9.11, 'C': 8.18, 'D': 9.60, 'E': 9.67, 'L': 8.95,
			 'R': 9.04, 'H': 9.17 }

	pKa = { 'D': 3.9, 'E': 4.3, 'R': 12.0, 'K': 10.5, 'H': 6.08, 'C': 8.28, 
			'Y': 10.1 }

	def pI( self, pH ):
		

