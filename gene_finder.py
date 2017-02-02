# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Camille Xue

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_comp = ''
    for i in reversed(dna):
        reverse_comp = reverse_comp + get_complement(i)
    return reverse_comp


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop_codons = ['TGA', 'TAA', 'TAG']
    for i in range(0, len(dna), 3):  # goes in steps of three through dna
        codon = dna[i:i+3]  # codon is set of three
        if codon in stop_codons:  # checks if it's in the stop_codons list
            return dna[0:i]
        else:
            continue
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    all_orfs = []
    end = -1  # for keeping track of the end of an ORF
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if codon == 'ATG' and i > end:  # makes sure non-nested by making sure it's after the end of the last orf
            all_orfs.append(rest_of_ORF(dna[i:len(dna)]))  # adds the orf to the list
            end = len(str(rest_of_ORF(dna[i:len(dna)])))  # marks end of orf
        else:
            continue
    return all_orfs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_orfs_all_frames = []
    all_orfs_all_frames.extend(find_all_ORFs_oneframe(dna))
    all_orfs_all_frames.extend(find_all_ORFs_oneframe(dna[1:]))
    all_orfs_all_frames.extend(find_all_ORFs_oneframe(dna[2:]))

    return all_orfs_all_frames
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_orfs_both_strands = []
    all_orfs_both_strands.extend(find_all_ORFs(dna))
    all_orfs_both_strands.extend(find_all_ORFs(get_reverse_complement(dna)))

    return all_orfs_both_strands
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_orf = find_all_ORFs_both_strands(dna)
    long_orf = max(all_orf, key=len)  # finds longest string in all_orf
    return long_orf
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    all_longest = []
    for n in range(0, num_trials):
        dna = shuffle_string(dna)
        all_longest.append(longest_ORF(dna))
    return len(max(all_longest))

    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acid = []
    if len(dna) % 3 == 1:
        dna = dna[:-1]
    elif len(dna) % 3 == 2:
        dna = dna[0:-2]

    for i in range(0, len(dna), 3):  # goes in steps of three through dna
        codon = dna[i:i+3]  # codon is set of three
        amino_acid.append(aa_table[codon])
    return ''.join(amino_acid)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    aa_strands = []
    threshhold = longest_ORF_noncoding(dna, 1500)
    all_ORF = find_all_ORFs_both_strands(dna)
    for i in all_ORF:
        if len(i) >= threshhold:
            aa_strands.append(coding_strand_to_AA(i))
        else:
            continue
    return aa_strands


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

    #doctest.run_docstring_examples(get_complement, globals(), verbose=True)
    #doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=True)
    ##doctest.run_docstring_examples(rest_of_ORF, globals(), verbose=True)
    #doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)
    #doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)
    #doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)
    #doctest.run_docstring_examples(longest_ORF, globals(), verbose=True)
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
