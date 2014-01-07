from yahmm import *
import itertools as it
import string
import numpy as np

'''
This module provides HMM based sequence-sequence, sequence-profile aligners for both pairwise
alignment and multiple sequence alignment.
'''

class PSSM( object ):
	'''
	A position specific scoring matrix. For every position in a multiple sequence alignment,
	characterize the density of each posibility. Takes in a MSA in list form and returns a
	easy to use representation.
	'''

	def __init__( self, msa, pseudocount=1e-4 ):
		'''
		Upon receiving the MSA, transform it into the PSSM.
		'''

		if type(msa) is list and type(msa[0]) is not list:
			msa = [ msa ]

		self.pssm = []
		for position in it.izip( *msa ):
			fpos = filter( lambda x: x is not '-', position )
			if len(fpos) == 0:
				continue

			pos_counts = { char: fpos.count(char) for char in set(fpos) }
			pseudocounts = { key: pos_counts[key] + pseudocount if key in pos_counts else
								  pseudocount for key in string.ascii_uppercase }

			total = sum( pseudocounts.values() ) * 1.
			self.pssm.append( { key: val / total for key, val in pseudocounts.items() } )

	def __getitem__( self, slice ):
		'''
		When one is slicing the PSSM, they really want the underlying data in the PSSM.
		'''

		return self.pssm[ slice ]

	def __repr__( self ):
		'''
		A string representation of the PSSM
		'''

		return '\n'.join( '\t'.join( "{} : {:2.2}".format( key, value ) for key, value in pos ) for pos in self.pssm )

	def __len__( self ):
		'''
		The number of positions in the PSSM.
		'''

		return len( self.pssm )

	def consensus( self ):
		return [ { value: key for key, value in position.items() }[ max( position.values() )] for position in self ]


class ProfileAligner( object ):
	def __init__( self, master, slave ):
		'''
		Must take in a PSSM object or a list of whatever is being aligned. Both x and y must be a
		list of at least one list, where each inner list represents a sequence. For example:

		x = [ [ 'A', 'B', 'C', 'D', 'E' ],
              [ '-', '-', 'C', 'D', 'E' ],
              [ 'A', 'B', 'D', '-', '-' ] ]
        y = [ [ 'A', 'B', 'E', 'F', 'G' ] ]

        This means that simple pairwise comparison can be done by generating PSSMs where each
        character has a ~100% probability in its respective position. All alignments are
        generalized as profile alignments in this manner.
		'''

		self.master = PSSM( master )
		self.slave = PSSM( slave )

	def _build( self, pssm ):
		model = Model( name="Global Profile Aligner" )
		insert_dist = { char: 1. / 26 for char in string.ascii_uppercase }

		last_match = model.start
		last_insert = State( DiscreteDistribution( insert_dist ), name="I0" )
		last_delete = None

		model.add_transition( model.start, last_insert, 0.15 )
		model.add_transition( last_insert, last_insert, 0.20 )

		for i, position in enumerate( pssm ):
			match = State( DiscreteDistribution( position ), name="M"+str(i+1) ) 
			insert = State( DiscreteDistribution( insert_dist ), name="I"+str(i+1) )
			delete = State( None, name="D"+str(i+1) )

			model.add_transition( last_match, match, 0.60 )
			model.add_transition( last_match, delete, 0.25 )
			model.add_transition( last_insert, match, 0.60 )
			model.add_transition( last_insert, delete, 0.20 )
			model.add_transition( delete, insert, 0.15 )
			model.add_transition( insert, insert, 0.20 )
			model.add_transition( match, insert, 0.15 )

			if last_delete:
				model.add_transition( last_delete, match, 0.60 )
				model.add_transition( last_delete, delete, 0.25 )

			last_match, last_insert, last_delete = match, insert, delete

		model.add_transition( last_delete, model.end, 0.80 )
		model.add_transition( last_insert, model.end, 0.80 )
		model.add_transition( last_match, model.end, 0.85 )

		model.bake()
		return model

	def global_alignment( self ):
		'''
		Perform a global alignment using a HMM. This aligns two profiles two each other,
		returning the probability of the alignment, and the two consensus alignments. 
		'''

		profile = self._build( self.master )
		prob, states = profile.viterbi( self.slave.consensus() )
		
		master = self.master.consensus()
		slave = self.slave.consensus()

		slave_aligned = []
		master_aligned = []
		i, j = 0, 0 
		# Follow y's path through x, ignoring start and end state
		for state in states[1:-1]:
			sname = state[1].name

			if sname.startswith( 'M' ):
				slave_aligned.append( slave[i] )
				master_aligned.append( master[j] )
				i += 1
				j += 1
			elif sname.startswith( 'D' ):
				slave_aligned.append( '-' )
				master_aligned.append( master[j] )
				j += 1
			elif sname.startswith( 'I' ):
				slave_aligned.append( slave[i] )
				master_aligned.append( '-' )  
				i += 1

		return prob, master_aligned, slave_aligned

class MultipleSequenceAligner( object ):
	def __init__( self, sequences ):
		self.sequences = sequences

	def _score( self, msa ):
		n = len( msa ) * 1.

		log_prob = lambda x, n: -x * math.log( x / n, 2 )
		entropy = lambda col: sum( log_prob( col.count( char ), n ) for char in set( col ) )  

		return sum( entropy( col ) for col in it.izip( *msa ) )

	def iterative_alignment( self, epsilon=3 ):
		'''
		Perform a HMM-based iterative alignment. If an initial alignment is provided, will use that
		to begin with, otherwise will simply use the sequences provided raw. This method will peel
		the top sequence off and align it to a profile of the other sequences, continuing this method
		until there is little change in the score for a full round of iteration. The scoring mechanism
		is done by minimum entropy.
		'''

		# Unpack the initial sequences
		msa = self.sequences

		n = len( msa )
		# Give initial scores
		score, last_score = float("-inf"), float("inf")

		# Until the scores converge...
		while abs( score-last_score ) >= epsilon:
			last_score = score
			score = 0
			# Run a full round of popping from the msa queue and enqueueing at the end
			for i in xrange( n ):
				# Pull a single 'master' sequence off the top of the msa queue
				master = filter( lambda x: x is not '-', msa[0] )

				# Make the rest of them slaves
				slaves = msa[1:]

				# Perform the alignment using the HMM-based profile aligner
				p, x, y = ProfileAligner( master, slaves ).global_alignment()

				# Add gaps to all sequences in the slave profile
				slaves.append( x )
				msa = slaves

				m = max( map( len, msa ) )
				for seq in msa:
					seq.extend( ['-']*(n-len(seq) ) )
					print seq
				print

				# Calculate the score for this run
				score += self._score( msa ) / float(n)
			print score

		return score, msa

if __name__ == '__main__':
	sequence = list('ABCDEFCDEFCDEFGCDEFHGIJK')

	msa = NaiveMSA( sequence )

	aligner = MultipleSequenceAligner( msa )
	score, msa = aligner.iterative_alignment()

	print score
	for seq in msa:
		print seq
