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

import argparse
import sys

from crumbs.utils.file_utils import get_input_fhand

# Missing docstring
# pylint: disable=C0111


def setup_basic_argparse(**kwargs):
    'It prepares the command line argument parsing.'

    parser = argparse.ArgumentParser(**kwargs)
    in_help = 'Input VCF file (default STDIN)'
    parser.add_argument('input', help=in_help, nargs='?',
                        type=argparse.FileType('rb'), default=sys.stdin)
    parser.add_argument('-o', '--output', default=sys.stdout,
                        help='Output VCF file (default STDOUT)',
                        type=argparse.FileType('w'))
    msg = 'File to print some statistics (default STDERR)'
    parser.add_argument('-l', '--log', help=msg, type=argparse.FileType('w'),
                        default=sys.stderr)

    return parser


def setup_filter_argparse(**kwargs):
    'It prepares the command line argument parsing.'
    parser = setup_basic_argparse(**kwargs)
    parser.add_argument('-f', '--filtered',
                        help='Output for filtered SNVs',
                        type=argparse.FileType('w'))
    parser.add_argument('-s', '--samples', action='append',
                        help='samples to use')
    parser.add_argument('-p', '--samples_file',
                        help='File with samples to use. One per line',
                        type=argparse.FileType('r'))
    return parser


def parse_basic_args(parser):
    parsed_args = parser.parse_args()

    in_fhand = get_input_fhand(parsed_args.input)

    out_fhand = parsed_args.output
    log_fhand = parsed_args.log

    args = {'in_fhand': in_fhand, 'log_fhand': log_fhand,
            'out_fhand': out_fhand}

    return args, parsed_args


def parse_sample_file(fhand):
    samples = []
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        samples.append(line)
    return samples


def parse_filter_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    filter_snvs_args, parsed_args = parse_basic_args(parser)

    filtered_fhand = parsed_args.filtered
    filter_snvs_args['filtered_fhand'] = filtered_fhand

    samples = set()
    if parsed_args.samples is not None:
        samples.update(parsed_args.samples)
    if parsed_args.samples_file is not None:
        samples.update(parse_sample_file(parsed_args.samples_file))

    filter_class_kwargs = {'samples_to_consider': samples if samples else None}

    return filter_snvs_args, filter_class_kwargs, parsed_args
