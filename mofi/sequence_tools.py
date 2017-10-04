"""
Functions/classes to handle protein sequences.
"""

import re
from mofi import mass_tools
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
    The chain number equals the number of header lines (starting with ``>``).
    The sequences of all chains are merged.

    :param str fasta_string: String in FASTA format
    :return: (1) number of chains,
             (2) string containing the raw sequence
    :rtype: tuple(int, str)
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
    Calculates the atoms in an amino acid sequence.

    :param str sequence: String of amino acids in one-letter format
    :param int chains: number of chains; for each chain, the weight of a
                       water molecule must be added to the total weight
    :param int disulfide_bonds: number of disulfide bonds
    :return: the elemental composition of the chains
             (example: ``{"C": 100; "H": 50; "N": 20}``)
    :rtype: dict
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

    :param str sequence: String containing an amino acid sequence
    :return: Two lists of N- and O-glycosylation sites, respectively.
             Each site is a tuple (site, (start, end)).
    :rtype: tuple(list, list)
    """
    n_pattern = re.compile("[ST]|N[^P][ST]")
    matchlist = [(m.group(), m.span()) for m in n_pattern.finditer(sequence)]
    n_site_list = [i for i in matchlist if len(i[0]) > 1]
    o_site_list = [i for i in matchlist if len(i[0]) == 1]
    return n_site_list, o_site_list


def apply_pngasef(sequence):
    """
    Simulate a digest with PNGase F, that is, change all asparagines
    in N-glycosylation sites to aspartates.

    :param str sequence: String of amino acids in one-letter format
    :return: (1) Sequence after treatment with PNGase F,
             (2) found N-glycan sites as returned by
                 :func:`find_glycosylation_sites()`
    :rtype: str, list(int)
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

    :ivar dict amino_acid_composition: {amino acid: count}
    :ivar mass_tools.Formula formula: atomic formula of the protein
    :ivar int n_sites: number of N-glycosylation sites
    :ivar float mass: the protein's molecular mass
    :ivar mass_without_disulfides: the reduced protein's molecular mass

    .. automethod:: __init__
    """

    def __init__(self, sequence, chains=1, disulfides=0, pngasef=False):
        """
        Create a new protein instance.

        :param str sequence: sequence of the protein
        :param int chains: number of chains
        :param int disulfides: number of disulfide bonds
        :param bool pngasef: true if the protein was treated with PNGase F
        """

        digested_sequence, self.n_sites = apply_pngasef(sequence)
        if pngasef:
            sequence = digested_sequence

        self.amino_acid_composition = {a: sequence.count(a)
                                       for a in amino_acid_names}
        self.formula = mass_tools.Formula(
            get_sequence_atoms(sequence, chains, disulfides))
        self.mass = self.formula.mass
        self.mass_without_disulfides = mass_tools.Formula(
            get_sequence_atoms(sequence, chains, 0)).mass
