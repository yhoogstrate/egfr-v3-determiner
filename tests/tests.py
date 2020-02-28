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


import unittest
import os
import pysam


def sam_to_sorted_bam(sam, sorted_bam):
    if os.path.exists(sorted_bam):
        return False
    else:
        # for save_stdout the file needs to be touched first
        fh = open(sorted_bam, 'wb')
        fh.close()
        
        # pysam in interactive mode requires the 'save_stdout' whereas in cli mode it requires '-o'
        pysam.sort('-o',sorted_bam, sam, save_stdout=sorted_bam)
        
        if os.path.getsize(sorted_bam) == 0:
            os.remove(sorted_bam)
            return False
        else:
            pysam.index(sorted_bam)

        return sorted_bam


TEST_DIR = "tests/data/"
TMP_DIR = "tmp/"


class Tests(unittest.TestCase):
    def test_001(self):
        input_file_sam = TEST_DIR + "test_001.sam"
        input_file_bam = TMP_DIR + "test_001.bam"
        
        sam_to_sorted_bam(input_file_sam, input_file_bam)
        
        from egfrviiideterminer import egfrviiideterminer
        dbkey = 'hg38'
        
        self.assertEqual(egfrviiideterminer.extract_viii_reads(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False), {'vIII': {'example_01'}, 'wt': set()})
        self.assertEqual(egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False), {'vIII': {'example_01'}, 'wt': set()})

    def test_002(self):
        input_file_sam = TEST_DIR + "test_002_wt_non-spliced.sam"
        input_file_bam = TMP_DIR + "test_002_wt_non-spliced.bam"
        
        sam_to_sorted_bam(input_file_sam, input_file_bam)
        
        from egfrviiideterminer import egfrviiideterminer
        dbkey = 'hg38'
        
        # hier moet die hem wel vinden, in wt
        self.assertEqual(egfrviiideterminer.extract_viii_reads(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False), {'vIII': set(), 'wt': {'example_002'}})
        
        # spliced only - hier moet die hem niet vinden, in wt
        self.assertEqual(egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False), {'vIII': set(), 'wt': {'example_002'}})

    def test_003(self):
        input_file_sam = TEST_DIR + "test_003_vIII_non_spliced.sam"
        input_file_bam = TMP_DIR + "test_003_vIII_non_spliced.bam"
        
        sam_to_sorted_bam(input_file_sam, input_file_bam)
        
        from egfrviiideterminer import egfrviiideterminer
        dbkey = 'hg19'
        
        # hier moet die hem wel vinden, in wt
        results = egfrviiideterminer.extract_viii_reads(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, True)
        self.assertEqual(len(results['vIII']), 170)
        self.assertEqual(len(results['wt']), 0)
        
        # do not allow PCR/optical duplicates
        results = egfrviiideterminer.extract_viii_reads(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False)
        self.assertEqual(len(results['vIII']), 43)
        self.assertEqual(len(results['wt']), 0)
        
        # spliced only - hier moet die hem niet vinden, in 
        results = egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, True)
        self.assertEqual(len(results['vIII']), 0)
        self.assertEqual(len(results['wt']), 0)

    def test_004(self):
        input_file_sam = TEST_DIR + "test_004_interchromosomal.sam"
        input_file_bam = TMP_DIR + "test_004_interchromosomal.bam"
        
        sam_to_sorted_bam(input_file_sam, input_file_bam)
        
        from egfrviiideterminer import egfrviiideterminer
        dbkey = 'hg19'
        
        # spliced only - hier moet die hem niet vinden, in 
        results = egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False)
        self.assertEqual(len(results['vIII']), 0)
        self.assertEqual(len(results['wt']), 0)

        # ook niet spliced, hier moet je hem wel in vinden
        results = egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], True, False)
        self.assertEqual(len(results['vIII']), 1)
        self.assertEqual(len(results['wt']), 0)

    def test_005(self):
        input_file_sam = TEST_DIR + "test_005_duplicate.sam"
        input_file_bam = TMP_DIR + "test_005_duplicate.bam"
        
        sam_to_sorted_bam(input_file_sam, input_file_bam)
        
        from egfrviiideterminer import egfrviiideterminer
        dbkey = 'hg19'
        
        # geen duplicates = niet vinden
        results = egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False)
        self.assertEqual(len(results['vIII']), 0)
        self.assertEqual(len(results['wt']), 1)

        # wel duplicates = wel vinden
        results = egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, True)
        self.assertEqual(len(results['vIII']), 0)
        self.assertEqual(len(results['wt']), 2)

    def test_006(self):
        input_file_sam = TEST_DIR + "test_006.sam"
        input_file_bam = TMP_DIR + "test_006.bam"
        
        sam_to_sorted_bam(input_file_sam, input_file_bam)
        
        from egfrviiideterminer import egfrviiideterminer
        dbkey = 'hg38'
        
        self.assertEqual(egfrviiideterminer.extract_viii_reads(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False), {'vIII': {'example_01'}, 'wt': set()})
        self.assertEqual(egfrviiideterminer.extract_viii_reads_based_on_sjs(input_file_bam, egfrviiideterminer.egfr_exons[dbkey], False, False), {'vIII': {'example_01'}, 'wt': set()})


if __name__ == '__main__':
    main()
