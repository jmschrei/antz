# io.py
# Contact: Jacob Schreiber
#          jmschr@cs.washington.edu

'''
This script focuses on data input and output, and currently supports the
following files:

* FastA
'''

from seq import *

class FastA( object ):
	'''
	This is a FastA file. It can contain many DNA, RNA, or Protein
	sequences in it. This can be read in or written out.
	'''

	def __init__( self, sequences ):
		'''
		If sequences are passed in, they should be as the DNA, RNA, or protein
		objects, so that all metadata is written out as well.
		'''

		self.sequences = sequences

	def __str__( self ):
		'''
		String representation of the FastA
		'''

		return '\n'.join( sequence.to_fasta() for sequence in self.sequences )


	def to_file( self, filename, attrs=None ):
		'''
		Write out a FastA file. Attrs specifies the attributes you want to
		write out as well, in that order. Since any data can be stored in these
		objects, it allows you to pick both what you want to write out, and
		in what order. If nothing is provided, nothing is written out.
		'''

		with open( filename, 'w' ) as outfile: 
			# Write out each stored sequence
			for sequence in self.sequences:
				outfile.write( sequence.to_fasta( attrs ) )

	@classmethod
	def from_file( cls, filename, attrs=None, delimiter=' ', seqType=None ):
		'''
		Read in a FastA file. Given names for each delimited item in the
		comments by specifying their attribute in order. Specify the seqType
		as the class object or string.
		'''

		if isinstance( seqType, str ):
			if seqType.lower() == 'protein':
				seqType = Protein
			elif seqType.lower() == 'rna':
				seqType = RNA
			elif seqType.lower() == 'dna':
				seqType = DNA
			else:
				seqType = Sequence
		seqType = seqType or Sequence

		sequences = []
		with open( filename, 'r' ) as infile:
			comments, sequence = None, ''

			# Go through the file line by line
			for line in infile:
				# If the next line starts with a >, it means that the previous
				# sequence has come to an end.
				if line.startswith( '>' ):
					# If a sequence has been found, create and append the
					# sequence object 
					if sequence != '':
						comments = comments.split( delimiter )
						attributes = { attr: comment for attr, comment in zip( attrs, comments ) }
						sequences.append( seqType( sequence, **attributes ) )
					
					# Now get the comment, removing the > and any newlines
					comments = line[1:].strip('\r\n')

					# Reset the sequence
					sequence = ''
				else:
					# Otherwise, append the sequence line to the growing
					# sequence
					sequence += line.strip('\r\n')

		comments = comments.split( delimiter )
		attributes = { attr: comment for attr, comment in zip( attrs, comments )}
		sequences.append( seqType( sequence, **attributes ) )
		return cls( sequences )