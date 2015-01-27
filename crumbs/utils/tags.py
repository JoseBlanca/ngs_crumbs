# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of ngs_crumbs.
# ngs_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# ngs_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ngs_crumbs. If not, see <http://www.gnu.org/licenses/>.


STDIN = 'stdin'
STDOUT = 'stdout'
INFILES = 'infiles'
OUTFILE = 'output'
# Input format tag when we want to guess
GUESS_FORMAT = 'guess'
PROCESSED_SEQS = 'processed_seqs'
PROCESSED_PACKETS = 'processed_packets'
YIELDED_SEQS = 'yielded_seqs'
UPPERCASE = 'upper'
LOWERCASE = 'lower'
SWAPCASE = 'swap'
FWD = 'fwd'
REV = 'rev'
TRIMMING_RECOMMENDATIONS = 'trimming_recommendations'
VECTOR = 'vector'
QUALITY = 'quality'
OTHER = 'other'
TRIMMING_KINDS = [VECTOR, QUALITY, OTHER]
ELONGATED = 'elongated'
SUBJECT = 'subject'
QUERY = 'query'
START = 0
END = 1
BGZF = 'bgzf'
GZIP = 'gzip'
BZIP2 = 'bzip2'

NUCL = 'nucl'
PROT = 'prot'

ERROR_ENVIRON_VARIABLE = 'ngs_crumbs_error_debugging'

FIVE_PRIME = "5'"
THREE_PRIME = "3'"

SEQS_PASSED = 'seqs_passed'
SEQS_FILTERED_OUT = 'seqs_filtered_out'

ORPHAN_SEQS = 'orphan_seqs'

SEQITEM = 'seqitem'
SEQRECORD = 'seqrecord'
SANGER_QUALITY = 'fastq'
ILLUMINA_QUALITY = 'fastq-illumina'
SANGER_FASTQ_FORMATS = ('fastq-sanger', 'fastq')
ILLUMINA_FASTQ_FORMATS = ('fastq-illumina',)

CHIMERA = 'chimera'
NON_CHIMERIC = 'non_chimeric'
UNKNOWN = 'unknown'
