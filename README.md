maptools
========

This package contains lift-over, alignment, and sequence retrieval tools

To install, type './configure' to cofigure and 'make all' to compile all the binaries

========
Utilities:

1. transf
    This utility transforms fasta files (usually, genomes) into an indexed 4-bit format,
    with 2 bits encoding ACGT and the other two encoding Ns and repeat masked state 
	
    ./transf -dir <dirname> -maskdir <dirname> -dbx <file> -idx <file> [-remove] [-uncompress] [-quiet]
	-dir name of the directory containing FASTA files
	-maskdir name of the directory containing masked FASTA files (optional)
	-maskext masked FASTA file extension
	-dbx, -idx sequence database names
	-remove remove source FASTA files [default=NO]
	-uncompress [default=NO]
	-quiet suppress verbose output [default=NO]

    Examples:
    ./transf -dir human_genome/file.fa -dbx hg19.dbx -idx hg19.idx
	Creates hg19.dbx and hg19.idx from all files in human_genome/
    ./transf -dir human_genome/file.fa -dbx hg19.dbx -idx hg19.idx -uncompress
	Creates human_genome/* from hg19.dbx and hg19.idx

2. getsegm
    This utility does sequence retrieval from the 4-bit repository (see transf) given the input file of intervals

    ./getsegm -in <file> -dbx <file> -idx <file> -out <file> [-limit <length>] [-margins <margin1> <margin2>] [-quiet] 
	-in BED file (NOTE: although the columns are as in BED file, the coordinate system is 1-based!)
 	-dbx, -idx database files
	-out output file
	-limit [default=800000] sequence length limit (if the length of the interval is greater than this limit, 
	nucleotides in the middle are replaced by the spacer)
	-spacer [default=.....]
	-margins <margin1> <margin2> [default=0 0] margins to add on each side of the interval

3. getwind
	This utility does sequence retrieval from the 4-bit repository (see transf) given the input file of positions and window sizes

    ./getwind -in <file> -dbx <file> -idx <file> -out <file>
	-we exonic window [default=0] -wi intronic window [default=150] -cis [use colunms 1-3] [default=4207174]
	-quiet <suppress verbose output> [default=no]
	-all <include all sites>
	-coord <offset for acceptor sites> [default=0]

