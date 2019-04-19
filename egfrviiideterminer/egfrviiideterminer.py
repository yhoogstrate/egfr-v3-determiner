#!/usr/bin/env python
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

"""[License: GNU General Public License v3 (GPLv3)]

    EGFR vIII determiner: counts vIII / non-vIII spliced reads in BAM files
    Copyright (C) 2019  Youri Hoogstrate

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    You can contact me via the github repository at the following url:
    <https://github.com/yhoogstrate/egfr-v3-determiner>

    You can e-mail me via 'y.hoogstrate' at the following webmail domain:
    gmail dot com
"""

import sys
import pysam

"""
Possibilities:

non-split:
[             ]       [             ]
      [====>--------------<====]


split-A:
[             ]       [             ]
           [== ....... ==>--<====]


split-B:
[             ]       [             ]
  [====>---<== ....... ==]


split-C:
[             ]       [             ]
           <== ....... ==]
           [== ....... ==>

EGVRvIII = splice variant van exon1 naar exon8

door insert size kan dit natuurlijk ook een read van exon1 naar exon9 of misschien zelfs exon10

max insert size moet dus uiteindelijk een parameter worden


rekenkundige truc:
vindt alle exon 1 reads die:
 - een deel alignen naar exons 2,3,4,5,6,7
 - een deel alignen naar exons 8,9,10 etc
"""

# hg19
egfr_exons = {
  'hg19': {
    "1": ['chr7',55086714,55087058],
    "2": ['chr7',55209979,55210130],
    "3": ['chr7',55210998,55211181],
    "4": ['chr7',55214299,55214433],
    "5": ['chr7',55218987,55219055],
    "6": ['chr7',55220239,55220357],
    "7": ['chr7',55221704,55221845],
    "8": ['chr7',55223523,55223639],
    "9": ['chr7',55224226,55224352],
    "10": ['chr7',55224452,55224525]
    },
  'hg38': {
    "1": ['chr7',55019021,55019365],
    "2": ['chr7',55142286,55142437],
    "3": ['chr7',55143305,55143488],
    "4": ['chr7',55146606,55146740],
    "5": ['chr7',55151294,55151362],
    "6": ['chr7',55152546,55152664],
    "7": ['chr7',55154011,55154152],
    "8": ['chr7',55155830,55155946],
    "9": ['chr7',55156533,55156659],
    "10": ['chr7',55156759,55156832]
    }
}



def get_splice_junction_positions(alignedsegment):
    """
    https://sourceforge.net/p/samtools/mailman/message/29373646/

    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    out = []
    #out = [None, None]# start, end
    offset = alignedsegment.reference_start
    for cigar in alignedsegment.cigartuples:
        s = [offset,cigar]
        
        # Insertion adds a gap, does not pad
        # soft and hardclip don't add padding
        #               M  D  =  X
        if cigar[0] in [0, 2, 7, 8]:# Match, Insertion, padding, Read match?, mismatch
            offset += cigar[1]
            s.append(offset)
        elif cigar[0] == 9:
            offset -= cigar[1]
        elif cigar[0] == 3: # splice junction starts
            #out[0] = offset
            #out[1] = offset + cigar[1] + 1

            out.append([offset, offset + cigar[1] + 1])
            
            offset += cigar[1]
    
    return out



def extract_viii_reads(bam, exons):
    set_2_7 = set([])
    set_8_10 = set([])
    readnames = {'1': set(),
                 '2': set_2_7, '3': set_2_7, '4': set_2_7, '5': set_2_7, '6': set_2_7, '7': set_2_7,
                 '8': set_8_10, '9': set_8_10, '10': set_8_10}

    fh = pysam.AlignmentFile(bam, "rb")

    for exon in exons:
        for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
            if read.get_overlap(exons[exon][1], exons[exon][2]):
                readnames[exon].add(read.query_name)

    total_intersection = readnames['1'].intersection(set_2_7, set_8_10)
    #if len(total_intersection) > 0:
    for _ in total_intersection:
        print( "Warning, read found aligned to exon1, one of the exons 2-7 AND one of the exons 8-10: " + _, file=sys.stderr)
    readnames['1'] = readnames['1'].difference(total_intersection) # important step, imagine a read that is aligned to exon1, one of the exons 2-7 AND one of the exons 8-10, that needs to be excluded
    
    exon1_to_exon2_7 = readnames['1'].intersection(set_2_7)
    exon1_to_exon8_10 = readnames['1'].intersection(set_8_10)

    return {'vIII': exon1_to_exon8_10, 'wt': exon1_to_exon2_7}


def extract_viii_reads_based_on_sjs(bam, exons):
    set_2 = set()# readnames of those that splice from exon 1 to 2
    set_8 = set()# readnames of those that splice from exon 1 to 8
    
    read_idx = {'wt': set_2, 'vIII': set_8}
    for exon in ['2', '8']:
            if int(exon) < 8:
                read_idx[exons[exon][1]] = set_2
            else:
                read_idx[exons[exon][1]] = set_8

    fh = pysam.AlignmentFile(bam, "rb")

    exon = "1"
    for read in fh.fetch(exons[exon][0], exons[exon][1], exons[exon][2]):
        if read.get_overlap(exons[exon][1], exons[exon][2]):
            for sj in get_splice_junction_positions(read):
                if sj[0] != None:
                    if sj[0] == exons['1'][2] and sj[1] in read_idx:
                        read_idx[sj[1]].add(read.query_name)
                        break

    #print(read_idx['vIII'])
    #print(read_idx['wt'])

    total_intersection = set_2.intersection(set_8)

    for _ in total_intersection:
        print( "Warning, read found aligned to exon1, one of the exons 2-7 AND one of the exons 8-10: " + _, file=sys.stderr)
    
    return {'vIII': read_idx['vIII'], 'wt': read_idx['wt']}
