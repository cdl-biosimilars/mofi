"""
    sequence_tools.py

    Helper functions/classes for protein sequences.

    Authors: Stefan Senn, Wolfgang Skala

    (c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

import re
import mass_tools
from collections import Counter


amino_acid_names = {
    "A": ("Alanine", "Ala"),
    "C": ("Cysteine", "Cys"),
    "D": ("Aspartic Acid", "Asp"),
    "E": ("Glutamic Acid", "Glu"),
    "F": ("Phenylalanine", "Phe"),
    "G": ("Glycine", "Gly"),
    "H": ("Histidine", "His"),
    "I": ("Isoleucine", "Ile"),
    "K": ("Lysine", "Lys"),
    "L": ("Leucine", "Leu"),
    "M": ("Methionine", "Met"),
    "N": ("Asparagine", "Asn"),
    "P": ("Proline", "Pro"),
    "Q": ("Glutamine", "Gln"),
    "R": ("Arginine", "Arg"),
    "S": ("Serine", "Ser"),
    "T": ("Threonine", "Thr"),
    "V": ("Valine", "Val"),
    "W": ("Tryptophane", "Trp"),
    "Y": ("Tyrosine", "Tyr"),
    "6": ("Pyroglutamic Acid", "Pyr"), 
    "7": ("Hydroxyproline", "Hyp"), 
    "J": ("Phosphoserine", "Sep"), 
    "Z": ("Phosphothreonine", "Tpo")}  

amino_acid_compositions = {
    "A": {"H": 5, "C": 3, "O": 1, "N": 1},
    "C": {"H": 5, "C": 3, "S": 1, "O": 1, "N": 1},
    "D": {"H": 5, "C": 4, "O": 3, "N": 1},
    "E": {"H": 7, "C": 5, "O": 3, "N": 1},
    "F": {"H": 9, "C": 9, "O": 1, "N": 1},
    "G": {"H": 3, "C": 2, "O": 1, "N": 1},
    "H": {"H": 7, "C": 6, "N": 3, "O": 1},
    "I": {"H": 11, "C": 6, "O": 1, "N": 1},
    "K": {"H": 12, "C": 6, "N": 2, "O": 1},
    "L": {"H": 11, "C": 6, "O": 1, "N": 1},
    "M": {"H": 9, "C": 5, "S": 1, "O": 1, "N": 1},
    "N": {"H": 6, "C": 4, "O": 2, "N": 2},
    "P": {"H": 7, "C": 5, "O": 1, "N": 1},
    "Q": {"H": 8, "C": 5, "O": 2, "N": 2},
    "R": {"H": 12, "C": 6, "N": 4, "O": 1},
    "S": {"H": 5, "C": 3, "O": 2, "N": 1},
    "T": {"H": 7, "C": 4, "O": 2, "N": 1},
    "V": {"H": 9, "C": 5, "O": 1, "N": 1},
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},
    "Y": {"H": 9, "C": 9, "O": 2, "N": 1},
    "6": {"H": 5, "C": 5, "O": 2, "N": 1},          # pyroglutamate
    "7": {"H": 7, "C": 5, "O": 2, "N": 1},          # hydroxyproline
    "J": {"H": 6, "C": 3, "O": 5, "N": 1, "P": 1},  # phosphoserine
    "Z": {"H": 8, "C": 4, "O": 5, "N": 1, "P": 1}}  # phosphothreonine


def read_fasta_string(fasta_string):
    """
    Extracts the chain number and sequence from a string in FASTA format.
    The chain number equals the number of header lines (i.e., lines starting with "<").
    The sequences of all chains are merged
    :param fasta_string: String in FASTA format
    :return: (1) number of chains, (2) string containing the raw sequence
    """
    seqences = []
    chains = 0
    for line in fasta_string.split("\n"):
        if line.startswith(">"):
            chains += 1
        else:
            seqences.append(line.strip())
    return chains, "".join(seqences)


def get_sequence_atoms(sequence, chains=1, disulfide_bonds=0):
    """
    Calculates the atoms in an amino acid sequence
    :param sequence: String of amino acids in one-letter format
    :param chains: number of chains; for each chain, the weight of a water molecule must be added to the total weight.
    :param disulfide_bonds: number of disulfide bonds
    :return: the elemental composition of the chains, represented as a dictionary
             (example: {"C": 100; "H": 50; "N": 20})
    """
    composition = {a: 0 for a in "CHNOPS"}
    sequence_aa_composition = Counter(sequence)
    for aa, count in sequence_aa_composition.items():
        for atom, atom_count in amino_acid_compositions[aa].items():
            composition[atom] += count * atom_count
    composition["H"] += chains * 2
    composition["O"] += chains * 1
    composition["H"] -= disulfide_bonds * 2
    return composition


def find_glycosylation_sites(sequence):
    """
    Searches a sequence for N or O-linked glycosylation sites.
    :param sequence: String containing an amino acid sequence
    :return: Two lists of N- and O-glycosylation sites, respectively.
             Each site is a tuple (site, (start, end)).
    """
    n_pattern = re.compile("[ST]|N[^P][ST]")
    matchlist = [(m.group(), m.span()) for m in n_pattern.finditer(sequence)]
    n_site_list = [i for i in matchlist if len(i[0]) > 1]
    o_site_list = [i for i in matchlist if len(i[0]) == 1]
    return n_site_list, o_site_list


def apply_pngasef(sequence):
    """
    Simulate a digest with PNGase F, that is, change all asparagines in N-glycosylation sites to aspartates.
    :param sequence: String of amino acids in one-letter format
    :return: (1) Sequence after treatment with PNGase F,
             (2) found N-glycosylation sites as returned by find_glycosylation_sites
    """
    site_list = find_glycosylation_sites(sequence)[0]
    if site_list:
        seq_list = list(sequence)
        for site in site_list:
            seq_list[site[1][0]] = "D"
        return "".join(seq_list), site_list
    else:
        return sequence, []


class Protein:
    """
        A class which represents proteins with a single chain.

        Members:
            n_sites
            average_mass
            monoisotopic_mass
            amino_acid_composition
            formula
    """

    def __init__(self, sequence, chains=1, disulfides=0, pngasef=False):
        """
        Returns a new Protein
        :param sequence: sequence of the protein
        :param chains: number of chains
        :param disulfides: number of disulfide bonds
        :param pngasef: bool, indicates whether the protein was treates with PNGase F
        """
        if pngasef:
            sequence, self.n_sites = apply_pngasef(sequence)
        else:
            self.n_sites = find_glycosylation_sites(sequence)[0]
        self.formula = mass_tools.Formula(get_sequence_atoms(sequence, chains, disulfides))
        self.average_mass = self.formula.average_mass
        self.monoisotopic_mass = self.formula.monoisotopic_mass
        self.amino_acid_composition = {}
        for a in amino_acid_names:
            self.amino_acid_composition[a] = sequence.count(a)
