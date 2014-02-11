Copyright 2012 Dmitri Pervouchine (dp@crg.eu), Lab Roderic Guigo
Bioinformatics and Genomics Group @ Centre for Genomic Regulation
Parc de Recerca Biomedica: Dr. Aiguader, 88, 08003 Barcelona 

This file is a part of the 'maptools' package.
'maptools' package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

'maptools' package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with 'maptools' package.  If not, see <http://www.gnu.org/licenses/>.

============================================================================

DESCRIPTION

maptools package contains a number of lift-over, alignment, and sequence retrieval tools.

============================================================================

INSTALLATION

To install, type './configure' to cofigure and 'make all' to compile all the binaries

============================================================================
UTILITIES

* Sequence compressing/decompressing/retrieval

transf: this utility transforms fasta files (usually, genomes) into an indexed 4-bit format,
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

getsegm: this utility does sequence retrieval from the 4-bit repository (see transf) given the input file of intervals

	./getsegm -in <file> -dbx <file> -idx <file> -out <file> [-limit <length>] [-margins <margin1> <margin2>] [-quiet] 
	-in BED file (NOTE: although the columns are as in BED file, the coordinate system is 1-based!)
	 -dbx, -idx database files
	-out output file
	-limit [default=800000] sequence length limit (if the length of the interval is greater than this limit, 
	nucleotides in the middle are replaced by the spacer)
	-spacer [default=.....]
	-margins <margin1> <margin2> [default=0 0] margins to add on each side of the interval

getwind: this utility does sequence retrieval from the 4-bit repository (see transf) given the input file of positions and window sizes

	./getwind -in <file> -dbx <file> -idx <file> -out <file>
	-we exonic window [default=0] -wi intronic window [default=150] -cis [use colunms 1-3] [default=4207174]
	-quiet <suppress verbose output> [default=no]
	-all <include all sites>
	-coord <offset for acceptor sites> [default=0]

* Lift-over and mapping utilities

map_agnostic: this utility does liftOver of nucleotide coordinates (cps format) by  using chain alignment.

	./map_agnostic -in <cps_file> -chain <chain_alignment_file> -out <output_file> [-margin <margin>] [-quiet]
	-in chromosome/position/strand tab-delimited file (strand is +/-), has to be sorted by position in ascending order
	-chain chain_alignment_file (see UCSC genome browser)
	-out <output_file> [default=stdout]
	-margin margin length [default=0]
	-quiet suppress verbose output [default=NO]

	Input format: chromosome1/position1/strand1; chain species 1=>2
	Output format: chromosome1/position1/strand1/chromosome2/position2/strand2/chain_id

map_single: same as map_agnostic, but also takes into account gene information specified in column 4 of cps

	Input format: chromosome1/position1/strand1; chain species 1=>2
	Output format: chromosome1/position1/strand1/chromosome2/position2/strand2/gene/site/type/score

syntenic_filter: this utility takes a non-unique mapping in cps3+cps3 format and does ad-hoc filtering of the projected 
	coordinates by maximum synteny genome-wide

	./syntenic_filter -in <aln_file> -out <output_file> [-maxdepth <int>] [-threshold <double>] [-lendiff <diff>] [-quiet]
	-in cps3+cps3 file, remember to sort by position in ascending order
	-out <output_file> [default=stdout]
	-maxdepth <integer> how many preceding positions can be skipped [default=4]
	-threshold <double> max change of segment length, in percent [default=1.50]
	-lendiff <integer>, [default=100000]
	-quiet suppress verbose output [default=NO]
	Note: the mapping [x,x+dx] -> [y,y+dy] is OK if |dy-dx|/dx<threshold OR |dy-dx|<dlimit

best_match: same as syntenic_filter but takes cps3+cps6 output of map_single and looks for maximum synteny PER GENE

net_filter: currently deprecated but still retained in the package.

================================================================================================
