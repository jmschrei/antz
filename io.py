def readFastA( fasta_file ):
	if type( fasta_file ) == str:
		file = open( fasta_file, 'r' )
	elif type( fasta_file ) == file:
		file = fasta_file
	else:
		raise TypeError( "Must pass in either a file or a pointer to a file." )

	seq, comments = '', ''
	for line in file:
		if line.startswith( ">" ):
			if comments != '':
				yield FastA( seq, comments )
			seq, comments = '', line.strip( "\r\n")
		else:
			seq += line.strip( "\r\n" )
	yield FastA( seq, comments ) 

class FastA( object ):
	def __init__( self, sequence=None, comments=None ):
		self.sequence=str(sequence)
		self.comments=str(comments)

	def to_file( self, filename, append=False ):
		if append:
			file = open( filename, 'a' )
		else:
			file = open( filename, 'w' )
		file.write( ">" + self.comments + "\n" )
		for i in xrange( 0, len(self.sequence), 80 ):
			try:
				file.write( self.sequence[i:i+80] + "\n" )
			except IndexError:
				file.write( self.sequence[i:] + "\n" )
