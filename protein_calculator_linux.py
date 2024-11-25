#!/usr/bin/env python
# coding: utf-8
# usage: pro_phys.py

import argparse
import re

parser = argparse.ArgumentParser(description='This script is used for batch prediction of physicochemical properties of protein sequences.\n'
                                 'usage: python pro_phys.py -f seq.fasta -o result.txt')
parser.add_argument('-f', '--seq_file', help='Please input fasta file containing multiple sequences.', required=True)
parser.add_argument('-o', '--out_file', help='Please input file name of outfile(****.txt,no spaces).', required=True)
args = parser.parse_args()
file_name = args.seq_file
out_file = args.out_file

##Predefined dictionaries
protein_weights = {
    'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
    'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
    'M': 131.1926, 'N': 114.1038, 'P': 97.1167, 'Q': 128.1307, 'R': 156.1875,
    'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.1760}

monoisotopic_protein_weights = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
    'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
    'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
    'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}

amino_acid_formulas = {
    "A": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 0},
    "C": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 1},
    "D": {"C": 4, "H": 7, "N": 1, "O": 4, "S": 0},
    "E": {"C": 5, "H": 9, "N": 1, "O": 4, "S": 0},
    "F": {"C": 9, "H": 11, "N": 1, "O": 2, "S": 0},
    "G": {"C": 2, "H": 5, "N": 1, "O": 2, "S": 0},
    "H": {"C": 6, "H": 9, "N": 3, "O": 2, "S": 0},
    "I": {"C": 6, "H": 13, "N": 1, "O": 2, "S": 0},
    "K": {"C": 6, "H": 14, "N": 2, "O": 2, "S": 0},
    "L": {"C": 6, "H": 13, "N": 1, "O": 2, "S": 0},
    "M": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 1},
    "N": {"C": 4, "H": 8, "N": 2, "O": 3, "S": 0},
    "P": {"C": 5, "H": 9, "N": 1, "O": 2, "S": 0},
    "Q": {"C": 5, "H": 10, "N": 2, "O": 3, "S": 0},
    "R": {"C": 6, "H": 14, "N": 4, "O": 2, "S": 0},
    "S": {"C": 3, "H": 7, "N": 1, "O": 3, "S": 0},
    "T": {"C": 4, "H": 9, "N": 1, "O": 3, "S": 0},
    "V": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 0},
    "W": {"C": 11, "H": 12, "N": 2, "O": 2, "S": 0},
    "Y": {"C": 9, "H": 11, "N": 1, "O": 3, "S": 0}
}

kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
      "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
      "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

DIWV = {"A": {"A": 1.0, "C": 44.94, "E": 1.0, "D": -7.49,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": -7.49,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": 1.0, "Y": 1.0},
        "C": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 20.26,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 33.60,
              "K": 1.0, "M": 33.60, "L": 20.26, "N": 1.0,
              "Q": -6.54, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": 33.60, "W": 24.68, "V": -6.54, "Y": 1.0},
        "E": {"A": 1.0, "C": 44.94, "E": 33.60, "D": 20.26,
              "G": 1.0, "F": 1.0, "I": 20.26, "H": -6.54,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 20.26, "S": 20.26, "R": 1.0,
              "T": 1.0, "W": -14.03, "V": 1.0, "Y": 1.0},
        "D": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": 1.0, "F": -6.54, "I": 1.0, "H": 1.0,
              "K": -7.49, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 1.0, "S": 20.26, "R": -6.54,
              "T": -14.03, "W": 1.0, "V": 1.0, "Y": 1.0},
        "G": {"A": -7.49, "C": 1.0, "E": -6.54, "D": 1.0,
              "G": 13.34, "F": 1.0, "I": -7.49, "H": 1.0,
              "K": -7.49, "M": 1.0, "L": 1.0, "N": -7.49,
              "Q": 1.0, "P": 1.0, "S": 1.0, "R": 1.0,
              "T": -7.49, "W": 13.34, "V": 1.0, "Y": -7.49},
        "F": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 13.34,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": -14.03, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": 1.0, "Y": 33.601},
        "I": {"A": 1.0, "C": 1.0, "E": 44.94, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 13.34,
              "K": -7.49, "M": 1.0, "L": 20.26, "N": 1.0,
              "Q": 1.0, "P": -1.88, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": -7.49, "Y": 1.0},
        "H": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -9.37, "F": -9.37, "I": 44.94, "H": 1.0,
              "K": 24.68, "M": 1.0, "L": 1.0, "N": 24.68,
              "Q": 1.0, "P": -1.88, "S": 1.0, "R": 1.0,
              "T": -6.54, "W": -1.88, "V": 1.0, "Y": 44.94},
        "K": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -7.49, "F": 1.0, "I": -7.49, "H": 1.0,
              "K": 1.0, "M": 33.60, "L": -7.49, "N": 1.0,
              "Q": 24.64, "P": -6.54, "S": 1.0, "R": 33.60,
              "T": 1.0, "W": 1.0, "V": -7.49, "Y": 1.0},
        "M": {"A": 13.34, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 58.28,
              "K": 1.0, "M": -1.88, "L": 1.0, "N": 1.0,
              "Q": -6.54, "P": 44.94, "S": 44.94, "R": -6.54,
              "T": -1.88, "W": 1.0, "V": 1.0, "Y": 24.68},
        "L": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": -7.49, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 33.60, "P": 20.26, "S": 1.0, "R": 20.26,
              "T": 1.0, "W": 24.68, "V": 1.0, "Y": 1.0},
        "N": {"A": 1.0, "C": -1.88, "E": 1.0, "D": 1.0,
              "G": -14.03, "F": -14.03, "I": 44.94, "H": 1.0,
              "K": 24.68, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": -6.54, "P": -1.88, "S": 1.0, "R": 1.0,
              "T": -7.49, "W": -9.37, "V": 1.0, "Y": 1.0},
        "Q": {"A": 1.0, "C": -6.54, "E": 20.26, "D": 20.26,
              "G": 1.0, "F": -6.54, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 20.26, "S": 44.94, "R": 1.0,
              "T": 1.0, "W": 1.0, "V": -6.54, "Y": -6.54},
        "P": {"A": 20.26, "C": -6.54, "E": 18.38, "D": -6.54,
              "G": 1.0, "F": 20.26, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": -6.54, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 20.26, "S": 20.26, "R": -6.54,
              "T": 1.0, "W": -1.88, "V": 20.26, "Y": 1.0},
        "S": {"A": 1.0, "C": 33.60, "E": 20.26, "D": 1.0,
              "G": 1.0, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 20.26, "P": 44.94, "S": 20.26, "R": 20.26,
              "T": 1.0, "W": 1.0, "V": 1.0, "Y": 1.0},
        "R": {"A": 1.0, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -7.49, "F": 1.0, "I": 1.0, "H": 20.26,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": 13.34,
              "Q": 20.26, "P": 20.26, "S": 44.94, "R": 58.28,
              "T": 1.0, "W": 58.28, "V": 1.0, "Y": -6.54},
        "T": {"A": 1.0, "C": 1.0, "E": 20.26, "D": 1.0,
              "G": -7.49, "F": 13.34, "I": 1.0, "H": 1.0,
              "K": 1.0, "M": 1.0, "L": 1.0, "N": -14.03,
              "Q": -6.54, "P": 1.0, "S": 1.0, "R": 1.0,
              "T": 1.0, "W": -14.03, "V": 1.0, "Y": 1.0},
        "W": {"A": -14.03, "C": 1.0, "E": 1.0, "D": 1.0,
              "G": -9.37, "F": 1.0, "I": 1.0, "H": 24.68,
              "K": 1.0, "M": 24.68, "L": 13.34, "N": 13.34,
              "Q": 1.0, "P": 1.0, "S": 1.0, "R": 1.0,
              "T": -14.03, "W": 1.0, "V": -7.49, "Y": 1.0},
        "V": {"A": 1.0, "C": 1.0, "E": 1.0, "D": -14.03,
              "G": -7.49, "F": 1.0, "I": 1.0, "H": 1.0,
              "K": -1.88, "M": 1.0, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 20.26, "S": 1.0, "R": 1.0,
              "T": -7.49, "W": 1.0, "V": 1.0, "Y": -6.54},
        "Y": {"A": 24.68, "C": 1.0, "E": -6.54, "D": 24.68,
              "G": -7.49, "F": 1.0, "I": 1.0, "H": 13.34,
              "K": 1.0, "M": 44.94, "L": 1.0, "N": 1.0,
              "Q": 1.0, "P": 13.34, "S": 1.0, "R": -15.91,
              "T": -7.49, "W": -9.37, "V": 1.0, "Y": 13.34},
        }

positive_pKs = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
negative_pKs = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
pKcterminal = {"D": 4.55, "E": 4.75}
pKnterminal = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7, }
charged_aas = ("K", "R", "H", "D", "E", "C", "Y")


class file_read:
    def __init__(self, filename):
        self.file = filename
    def pro_seq_id(self):
        pros = open(self.file, 'r')
        pro_seq_id = {}
        pro_info = []
        for fasta in pros:
            pro_info.append(str(fasta.replace('\n', '').split('*')[0]))
        pro_id_indexs = [index for index, item in enumerate(pro_info) if '>' in item]
        for pro_id_index in pro_id_indexs:
            pro_id = str(pro_info[pro_id_index].split(">")[1])
            if pro_id_indexs.index(int(pro_id_index)) < len(pro_id_indexs) - 1:
                end = int(pro_id_indexs[pro_id_indexs.index(int(pro_id_index)) + 1])
            else:
                end = len(pro_info)
            pro_seq = ''.join(pro_info[int(pro_id_index) + 1:end])
            pro_seq_id[pro_id] = pro_seq
        return pro_seq_id


class ProteinAnalysis:
    def __init__(self, protein):
        self.protein = protein.upper()
        self.amino_acids = list(self.protein)
        self.length = len(self.protein)

    def number_of_amine_acid(self):
        count_amine_acid = len(self.protein)
        return count_amine_acid

    def molecular_weight(self, monoisotopic=False):
        if not self.protein:
            return None
        water = 18.01056 if monoisotopic else 18.01524
        weight_table = monoisotopic_protein_weights if monoisotopic else protein_weights
        weight = sum(weight_table[aa] for aa in self.amino_acids)

        weight += water

        return weight

    def gravy(self):
        if not self.protein:
            return None
        total_gravy = sum(kd[aa] for aa in self.amino_acids)
        return total_gravy / len(self.protein)

    def instability_index(self):
        index = DIWV
        score = 0.0

        for i in range(self.length - 1):
            this, next = self.protein[i: i + 2]
            dipeptide_value = index[this][next]
            score += dipeptide_value
        return (10.0 / self.length) * score

    def calculate_protein_formula(self):
        formulas = amino_acid_formulas
        total_formula = {"C": 0, "H": 0, "N": 0, "O": 0, "S": 0}

        for amino_acid in self.protein:
            if amino_acid in formulas:
                for element, count in formulas[amino_acid].items():
                    total_formula[element] += count

        num_peptide_bonds = len(self.protein) - 1
        total_formula["H"] -= num_peptide_bonds * 2
        total_formula["O"] -= num_peptide_bonds
        total_formula["C"] += 1  # N-terminal amino group
        total_formula["H"] += 1  # C-terminal carboxyl group
        total_formula["N"] += 1  # N-terminal amino group
        total_formula["O"] += 1  # C-terminal carboxyl group
        total_formula["S"] += 1  # N-terminal amino group
        final_formula = ''.join([f"{element}{total_formula[element] - 1}" for element in sorted(total_formula.keys()) if
                                 total_formula[element] > 0])
        pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
        elements = pattern.findall(final_formula)
        atom_count = {}
        total_atoms = 0
        for element, count in elements:
            if count == '':
                count = 1
            else:
                count = int(count)
            atom_count[element] = count
            total_atoms += count
        return final_formula, total_atoms

    def calculate_aliphatic_index(self):
        aliphatic_index = 100 * (sum(1 for aa in self.protein if aa == 'A') / len(self.protein) +
                                 2.9 * sum(1 for aa in self.protein if aa == 'V') / len(self.protein) +
                                 3.9 * (sum(1 for aa in self.protein if aa == 'I') +
                                        sum(1 for aa in self.protein if aa == 'L')) / len(self.protein))
        return aliphatic_index


class IsoelectricPoint:
    def __init__(self, protein, aa_content=None):
        """Initialize the class."""
        self.protein = protein.upper()
        if not aa_content:
            from Bio.SeqUtils.ProtParam import ProteinAnalysis as _PA

            aa_content = _PA(self.protein).count_amino_acids()
        self.charged_aas_content = self._select_charged(aa_content)

        self.pos_pKs, self.neg_pKs = self._update_pKs_tables()

    def _select_charged(self, aa_content):
        charged = {}
        for aa in charged_aas:
            charged[aa] = float(aa_content[aa])
        charged["Nterm"] = 1.0
        charged["Cterm"] = 1.0
        return charged

    def _update_pKs_tables(self):
        pos_pKs = positive_pKs.copy()
        neg_pKs = negative_pKs.copy()
        nterm, cterm = self.protein[0], self.protein[-1]
        if nterm in pKnterminal:
            pos_pKs["Nterm"] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs["Cterm"] = pKcterminal[cterm]
        return pos_pKs, neg_pKs

    def charge_at_pH(self, pH):
        positive_charge = 0.0
        for aa, pK in self.pos_pKs.items():
            partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
            positive_charge += self.charged_aas_content[aa] * partial_charge

        negative_charge = 0.0
        for aa, pK in self.neg_pKs.items():
            partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
            negative_charge += self.charged_aas_content[aa] * partial_charge

        return positive_charge - negative_charge

    def pi(self, pH=7.775, min_=4.05, max_=12):
        charge = self.charge_at_pH(pH)
        if max_ - min_ > 0.0001:
            if charge > 0.0:
                min_ = pH
            else:
                max_ = pH
            next_pH = (min_ + max_) / 2
            return self.pi(next_pH, min_, max_)
        return pH


with open(out_file, "w", encoding='utf-8') as f:
    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        'seq_id',
        'number_of_amine_acid',
        'molecular_weight',
        'theoretical_pi',
        'Formula',
        'Total_number_of_atoms',
        'instability_index',
        'aliphatic_index',
        'gravy'))
    i = 0
    pros = file_read(file_name)
    for pro in pros.pro_seq_id():
        pro_id = pro
        pro_seq = pros.pro_seq_id()[pro_id]
        prot_param = ProteinAnalysis(pro_seq)
        prot_pi = IsoelectricPoint(pro_seq)
        print("=" * 10, "seq", i + 1, "->", pro_id, "on the way", "=" * 10)
        number_of_amine_acid = len(pro_seq)
        molecular_weight = prot_param.molecular_weight()
        theoretical_pi = prot_pi.pi()
        Formula = prot_param.calculate_protein_formula()[0]
        Total_number_of_atoms = prot_param.calculate_protein_formula()[1]
        instability_index = prot_param.instability_index()
        aliphatic_index = prot_param.calculate_aliphatic_index()
        gravy = prot_param.gravy()
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            pro_id,
            number_of_amine_acid,
            round(molecular_weight, 2),
            round(theoretical_pi, 2),
            Formula,
            Total_number_of_atoms,
            round(instability_index, 2),
            round(aliphatic_index, 2),
            round(gravy, 3)))
        i += 1
