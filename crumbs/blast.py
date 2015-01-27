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

import os.path
import subprocess
import tempfile

from crumbs.utils.optional_modules import NCBIWWW
from crumbs.seq.seqio import seqio, guess_seq_type, write_seqs
from crumbs.utils.bin_utils import (check_process_finishes, popen,
                                    get_binary_path)
from crumbs.utils.tags import NUCL, PROT
from crumbs.seq.alignment_result import (filter_alignments, ELONGATED, QUERY,
                                         covered_segments_from_match_parts,
                                         elongate_match_parts_till_global,
                                         TabularBlastParser, BlastParser)
from crumbs.utils.file_utils import TemporaryDir
from crumbs.settings import get_setting


BLAST_FIELDS = {'query': 'qseqid', 'subject': 'sseqid', 'identity': 'pident',
                'aligment_length': 'length', 'mismatches': 'mismatch',
                'gap_open': 'gapopen', 'query_start': 'qstart',
                'query_end': 'qend', 'subject_start': 'sstart',
                'subject_end': 'send', 'expect': 'evalue', 'score': 'bitscore',
                'subject_length': 'slen', 'query_length': 'qlen'}

TASKS = {'blastn': ('megablast', 'blastn', 'blastn-short', 'dc-megablast',
                   'rmblastn'),
        'blastp': ('blastp', 'blastp-short', 'deltablast')}

REMOTE_BLAST_DBS = ['nt', 'nr']


def generate_tabblast_format(fmt):
    'Given a list with fields with our names it return one with the blast ones'
    return '6 ' + ' '.join([BLAST_FIELDS[f] for f in fmt])


def _makeblastdb_plus(seq_fpath, dbtype, outputdb=None):
    'It creates the blast db database'
    cmd = [get_binary_path('makeblastdb'), '-in', seq_fpath, '-dbtype', dbtype]
    if outputdb is not None:
        cmd.extend(['-out', outputdb])
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])


def _get_abs_blastdb_path(blastdb_or_path, dbtype):
    'It returns the blastdb absolute fpath'
    if os.path.isabs(blastdb_or_path):
        return blastdb_or_path
    else:
        paths = [blastdb_or_path]
        if 'BLASTDB' in os.environ:
            paths.append(os.path.join(os.environ['BLASTDB'], blastdb_or_path))
        for path in paths:
            if (os.path.exists(path) or _blastdb_exists(path, dbtype)):
                return os.path.abspath(path)
    return OSError('blastdb not found')


def _blastdb_exists(dbpath, dbtype=None):
    'It checks if a a blast db exists giving its path'
    assert dbtype in (NUCL, PROT, None)
    if dbtype == NUCL:
        exts = ['.nin', '.00.nin']
    elif dbtype == PROT:
        exts = ['.pin', '.00.pin']
    else:
        exts = ['.nin', '.pin', '.00.nin', '.00.pin']
    return any([os.path.exists(dbpath + ext) for ext in exts])


def get_or_create_blastdb(blastdb_or_path, dbtype=None, directory=None):
    '''it returns a blast database.

    If it does not exists it creates and if you give it a directory it will
    create it in that directory if it does not exist yet
    '''
    seq_fpath = _get_abs_blastdb_path(blastdb_or_path, dbtype)
    if directory:
        dbname = os.path.basename(seq_fpath)
        dbpath = os.path.join(directory, dbname)
    else:
        dbpath = seq_fpath

    if not _blastdb_exists(dbpath, dbtype=dbtype):
        if not os.path.exists(seq_fpath):
            msg = 'An input sequence is required to create a blastdb'
            raise RuntimeError(msg)
        if seq_fpath != dbpath:
            seqio([open(seq_fpath)], open(dbpath, 'w'), out_format='fasta',
                  copy_if_same_format=False)
        if dbtype is None:
            dbtype = guess_seq_type(open(dbpath))
        _makeblastdb_plus(dbpath, dbtype)
    return dbpath


def _do_blast_remote(query_fpath, database, program, out_fpath, params=None):
    'It does a blast using NCBI WWW server'
    # output always XML, although it can be configured the types are not
    # compatible with command line output types
    if not params:
        params = {}
    evalue, task = _parse_blast_params(params, program)
    if database not in REMOTE_BLAST_DBS:
        msg = "Can't perform a blast against {} in NCBI servers"
        msg = msg.format(database)
        raise RuntimeError(msg)
    if 'outfmt' in params and params['outfmt'] != "XML":
        msg = "To performs a remote blast the output format must be xml"
        raise RuntimeError(msg)

    ncbi_params = {'program': program, 'database': database,
                   'sequence': open(query_fpath).read(), 'expect': evalue}
    if task:
        ncbi_params['service'] = task

    result_handle = NCBIWWW.qblast(**ncbi_params)
    out_fhand = open(out_fpath, 'w')
    out_fhand.write(result_handle.read())
    out_fhand.close()


def _parse_blast_params(params, program):
    if 'expect' in params:
        if 'evalue' in params:
            msg = 'expect and evalue cannot be given as params at the same '
            msg += 'time'
            raise ValueError(msg)
        else:
            params['evalue'] = params['expect']
            del params['expect']
    if 'evalue' in params:
        evalue = params['evalue']
        del params['evalue']
    else:
        evalue = 0.001

    available_tasks = TASKS.get(program, None)
    if available_tasks:
        if 'task' in params:
            task = params['task']
            del params['task']
        else:
            task = available_tasks[0]
        assert task in available_tasks
    else:
        task = None

    return evalue, task


def do_blast(query_fpath, db_fpath, program, out_fpath, params=None,
             remote=False):
    if remote:
        _do_blast_remote(query_fpath, db_fpath, program, out_fpath,
                         params=params)
    else:
        _do_blast_local(query_fpath, db_fpath, program, out_fpath,
                        params=params)


def _do_blast_local(query_fpath, db_fpath, program, out_fpath, params=None):
    'It does a blast'
    if not params:
        params = {}
    evalue, task = _parse_blast_params(params, program)

    if 'outfmt' in params:
        outfmt = params['outfmt']
        del params['outfmt']
    else:
        outfmt = 5

    if program not in ('blastn', 'blastp', 'blastx', 'tblastx', 'tblastn'):
        raise ValueError('The given program is invalid: ' + str(program))
    binary = get_binary_path(program)
    cmd = [binary, '-query', query_fpath, '-db', db_fpath, '-out', out_fpath]
    cmd.extend(['-evalue', str(evalue), '-outfmt', str(outfmt)])
    if task:
        cmd.extend(['-task', task])
    if params:
        for key, value in params.viewitems():
            cmd.extend(('-' + key, str(value)))
    process = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    check_process_finishes(process, binary=cmd[0])


def _do_blast_2(db_fpath, queries, program, dbtype=None, blast_format=None,
                params=None, remote=False):
    '''It returns an alignment result with the blast.

    It is an alternative interface to the one based on fpaths.
    db_fpath should be a plain sequence file.
    queries should be a SeqRecord list.
    If an alternative blast output format is given it should be tabular, so
    blast_format is a list of fields.
    '''

    query_fhand = write_seqs(queries, file_format='fasta')
    query_fhand.flush()

    if remote:
        blastdb = db_fpath
        fmt = 'XML' if blast_format is None else blast_format.upper()
    else:
        blastdb = get_or_create_blastdb(db_fpath, dbtype=dbtype)
        if blast_format is None:
            blast_format = ['query', 'subject', 'query_length', 'subject_length',
                            'query_start', 'query_end', 'subject_start',
                            'subject_end', 'expect', 'identity', ]
        fmt = generate_tabblast_format(blast_format)

    if params is None:
        params = {}
    params['outfmt'] = fmt

    blast_fhand = tempfile.NamedTemporaryFile(suffix='.blast')
    do_blast(query_fhand.name, blastdb, program, blast_fhand.name, params,
             remote=remote)
    if remote:
        blasts = BlastParser(blast_fhand)
    else:
        blasts = TabularBlastParser(blast_fhand, blast_format)

    return blasts, blast_fhand


class BlasterForFewSubjects(object):
    '''It matches the given SeqRecords against the reads in the file.

    This class uses Blast to do the matching and it is optimized for having
    few SeqRecords to match against a medium sized sequence file.

    In this class the seqs_fpath, that in the API acts as the query, it's
    internally indexed into a blast database and the seqrecords are blasted
    against it. This is done because in the case of having lots of queries to
    blast against few subjects it is more efficient to do the blast backwards.
    This class behaves as if the seqs_fpath were the query and the seqrecords
    the subject although internally the blast is done with the query and
    subject changed.
    '''
    def __init__(self, seqs_fpath, seqs, program, params=None,
                 filters=None, elongate_for_global=False, seqs_type=None):
        '''It inits the class.'''
        self.program = program
        if params is None:
            params = {}
        params['max_target_seqs'] = str(get_setting('PACKET_SIZE'))
        self.params = params
        if filters is None:
            filters = []
        self.filters = filters
        self.elongate_for_global = elongate_for_global
        self._match_parts = self._look_for_blast_matches(seqs_fpath,
                                                         seqs,
                                                         seqs_type)

    def _look_for_blast_matches(self, seq_fpath, oligos, seqs_type):
        'It looks for the oligos in the given sequence files'
        # we need to keep the blast_fhands, because they're temp files and
        # otherwise they might be removed
        temp_dir = TemporaryDir()
        dbpath = os.path.join(temp_dir.name, os.path.basename(seq_fpath))
        seqio([open(seq_fpath)], open(dbpath, 'w'), out_format='fasta',
              copy_if_same_format=False)

        blasts, blast_fhand = _do_blast_2(dbpath, oligos, params=self.params,
                                          program=self.program,
                                          dbtype=seqs_type)
        if self.filters is not None:
            blasts = filter_alignments(blasts, config=self.filters)

        # Which are the regions covered in each sequence?
        indexed_match_parts = {}
        for blast in blasts:
            oligo = blast['query']
            for match in blast['matches']:
                read = match['subject']
                if self.elongate_for_global:
                    elongate_match_parts_till_global(match['match_parts'],
                                                 query_length=oligo['length'],
                                                 subject_length=read['length'],
                                                 align_completely=QUERY)

                # match_parts = [m['match_parts'] for m in blast['matches']]
                match_parts = match['match_parts']
                try:
                    indexed_match_parts[read['name']].extend(match_parts)
                except KeyError:
                    indexed_match_parts[read['name']] = match_parts

        temp_dir.close()
        blast_fhand.close()
        return indexed_match_parts

    def get_matched_segments_for_read(self, read_name):
        'It returns the matched segments for any oligo'
        setting_key = 'DEFAULT_IGNORE_ELONGATION_SHORTER'
        ignore_elongation_shorter = get_setting(setting_key)

        try:
            match_parts = self._match_parts[read_name]
        except KeyError:
            # There was no match in the blast
            return None

        # Any of the match_parts has been elongated?
        elongated_match = False
        for m_p in match_parts:
            if ELONGATED in m_p and m_p[ELONGATED] > ignore_elongation_shorter:
                elongated_match = True
        segments = covered_segments_from_match_parts(match_parts,
                                                     in_query=False)
        return segments, elongated_match


class Blaster(object):
    '''It matches the given SeqRecords against a blast database.
    It needs iterable with seqrecords and a blast dabatase
    '''

    def __init__(self, seqrecords, blastdb, program, dbtype=None, params=None,
                 filters=None, remote=False):
        self.program = program
        if params is None:
            params = {}
        self.params = params
        if filters is None:
            filters = []
        self.filters = filters
        self._remote = remote
        if dbtype not in (NUCL, PROT, None):
            raise ValueError('dbtype must be NUCL, PROT or None (we guess)')
        self._blasts = self._look_for_blast_matches(seqrecords, blastdb,
                                                    dbtype)

    def _look_for_blast_matches(self, seqrecords, blastdb, dbtype):
        'it makes the blast and filters the results'
        blasts, blast_fhand = _do_blast_2(blastdb, seqrecords, self.program,
                                          params=self.params, dbtype=dbtype,
                                          remote=self._remote)
        # print open(blast_fhand.name).read()
        if self.filters is not None:
            blasts = filter_alignments(blasts, config=self.filters)

        blasts = {blast['query']['name']: blast for blast in blasts}
        blast_fhand.close()
        return blasts

    def get_matched_segments(self, seqrecord_name):
        'It returns the matched segments for the given query (SeqRecord) name.'
        try:
            blast = self._blasts[seqrecord_name]
        except KeyError:
            # There was no match in the blast
            return None

        match_parts = [mp for m in blast['matches'] for mp in m['match_parts']]
        segments = covered_segments_from_match_parts(match_parts,
                                                     in_query=False)
        return segments

    @property
    def blasts(self):
        'It returns the blasts alignment results.'
        # we trust the clients not to mess with this dict or its nested
        # structure, to deepcopy would be a waste of resources in a
        # well-behaved world
        return self._blasts
