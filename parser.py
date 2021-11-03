# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 13:51:13 2021

@author: uha18
"""

#from math import sqrt
#def main():
#negacids = ["ASP", "GLU"]
#posacids = ["LYS", "ARG", "HIS"]
#pos_count = 0
#neg_count = 0
#hetero_count = 0
#sequence = []
#def interatom_dist(atom1, atom2):
#    return math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
path = input("Please provide the path to the .pdb file to analyze:")
with open(path, 'r') as pdbfile:
    #aminoacids {"ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "GLY". "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    element_dict = dict() #creates the dictionary of elements to occur in the given structure
    aminoacid_dict = dict() # creates the dictionary of aminoacid to occur in the given structure
    lig_dict = dict() # creates the dictionary for ligands
    negacids = ["ASP", "GLU"] # creates and populates the list of codes of aminoacids with negatively charged side chains at physiological pH
    posacids = ["LYS", "ARG", "HIS"] # creates and populates the list of codes of aminoacids with positively charged side chains at physiological pH
    pos_count = 0
    neg_count = 0
    hetero_count = 0 
    sequence = []
    
    for line in pdbfile:
        strcont = line.split() #splitting the line into a list of strings  
        if strcont[0] == "ATOM":
            atom_type = strcont[2]
            chain_id = strcont[4]
            element_type = strcont[11]
            if element_type in element_dict:
                element_dict[element_type] += 1
            else:
                element_dict[element_type] = 1
            if atom_type == "CA":
                aminoacid_code = strcont[3]
                if aminoacid_code in aminoacid_dict:
                    aminoacid_dict[aminoacid_code] += 1
                else:
                    aminoacid_dict[aminoacid_code] = 1
                if aminoacid_code in negacids:
                    neg_count += 1
                elif aminoacid_code in posacids:
                    pos_count +=1
                sequence.append(aminoacid_code)
                #calpha_coord = (strcont[6], strcont[7], strcont[8])
            elif strcont[0] == "HETATM" and strcont[3] != "HOH":
                pass               
    #print(sequence)
    print("Amino acid composition:")
    for keyaa in list(aminoacid_dict.keys()):
        print(keyaa, aminoacid_dict[keyaa], str(float('%.2f' %(100*aminoacid_dict[keyaa]/len(sequence))))+'%')
    print("Atomic composition:")
    for keyel in list(element_dict.keys()):
        print(keyel, element_dict[keyel])
    print("Total number of positively charged residues:", pos_count)
    print("Total number of negatively charged residues:", neg_count)
    #print ("Total number of ligands:", len(list(lig_set())))
        