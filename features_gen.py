#SSF

# Secondary Structure Folds

# Libraries
import pandas as pd
import numpy as np
!pip3 install wget
import wget
import re
!apt-get install dssp
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import pickle
from sklearn.preprocessing import StandardScaler, LabelEncoder

# Tasks
# Separate window sizes (3, 5, 7, 9, 11, 13)
window = 9

# dataset import and preprocessing
ds = pd.read_csv('newdataset.csv')
pdb = list(ds.iloc[:,1].values)
new_pdb = []
for i in pdb:
	i = i.replace('.pdb', '')
	new_pdb.append(i)

pdb = new_pdb
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

# Structures
# H,G,I: 1 
# T: 2 (T, S)
# S: 3
# B: 4
# E: 5
# - 6
# Exception: 7

ssf_list = []
p = PDBParser()
last_file = '/content/pdb/' + str(pdb[0]) + '.pdb'
last_dssp = dssp_dict_from_pdb_file(last_file)
for i in range(len(pdb)):
	try:
		pdb_id = str(pdb[i])
		print(pdb_id)
		try:
			file = '/content/pdb/' + pdb_id.lower() +'.pdb'
			if file == last_file:
				dssp = last_dssp
			else:
				last_file = file 
				dssp = dssp_dict_from_pdb_file(file)
				last_dssp = dssp
		except:
			file = '/content/pdb/' + pdb_id.upper() +'.pdb'
			if file == last_file:
				dssp = last_dssp
			else:
				last_file = file 
				dssp = dssp_dict_from_pdb_file(file)
				last_dssp = dssp
		dssp = dssp[0]
		ssf = []
		start = res[i] - window
		end = res[i] + window
		structure = ''
		for k, v in dssp:
			chain = k
			break
		for j in range(start-1, end):
			try:
				structure = dssp[chain, (' ', j, ' ')][1]
				if structure == 'H' or structure == 'G' or structure == 'I':
					ssf.append(1)
				elif structure == 'T' or structure == 'S':
					ssf.append(2)
				elif structure == 'B':
					ssf.append(3)
				elif structure == 'E':
					ssf.append(4)
				else:
					ssf.append(5)
			except:
				ssf.append(6)
		print(ssf, i, len(pdb))
	except:
		print("Error")
		ssf = np.zeros(window*2 + 1, dtype = int)
		print(ssf, i, len(pdb))
	ssf_list.append(ssf)

# Pickle
filename = 'NF1_' + str(window) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(ssf_list,outfile)
outfile.close()


#AASIS

# Amino Acids Signatures in the interaction shells

!pip3 install biopython
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np
import pickle

# create parser
parser = PDBParser()

# parameters
radius = 6

# read structure from file
ds = pd.read_csv('newdataset.csv')
pdb = list(ds.iloc[:,1].values)
new_pdb = []
for i in pdb:
	i = i.replace('.pdb', '')
	new_pdb.append(i)

pdb = new_pdb
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

# Iterate for the protein
inshell_proteins_5 = []
inshell_proteins_6 = []
inshell_proteins_7 = []
inshell_proteins_8 = []
for i in range(len(pdb)):
	print(i, len(pdb))
	inshell_proteins_single_5 = np.zeros(20, dtype = int) 
	inshell_proteins_single_6 = np.zeros(20, dtype = int) 
	inshell_proteins_single_7 = np.zeros(20, dtype = int) 
	inshell_proteins_single_8 = np.zeros(20, dtype = int) 
	try:
		try:
			path_ = '/content/pdb/' + str(pdb[i]).lower() + '.pdb'
			structure = parser.get_structure('PHA-L', path_)
		except:
			path_ = '/content/pdb/' + str(pdb[i]).upper() + '.pdb'
			structure = parser.get_structure('PHA-L', path_)
		model = structure[0]
		try:
			for chain in model:
				residue1 = chain[resid[i]] 
				for residue2 in chain:
					if residue1 != residue2:
						try:
							distance = residue1['CA'] - residue2['CA']
						except KeyError:
							continue
						if distance < 5: 
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_5[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_5[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_5[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_5[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_5[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_5[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_5[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_5[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_5[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_5[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_5[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_5[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_5[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_5[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_5[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_5[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_5[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_5[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_5[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_5[19] += 1
						if distance < 6:
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_6[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_6[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_6[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_6[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_6[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_6[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_6[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_6[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_6[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_6[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_6[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_6[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_6[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_6[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_6[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_6[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_6[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_6[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_6[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_6[19] += 1
						if distance < 7:
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_7[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_7[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_7[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_7[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_7[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_7[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_7[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_7[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_7[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_7[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_7[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_7[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_7[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_7[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_7[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_7[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_7[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_7[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_7[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_7[19] += 1
						if distance < 8:
							if residue2.get_resname() == 'ALA':
								inshell_proteins_single_8[0] += 1
							elif residue2.get_resname() == 'ARG':
								inshell_proteins_single_8[1] += 1
							elif residue2.get_resname() == 'ASN':
								inshell_proteins_single_8[2] += 1
							elif residue2.get_resname() == 'ASP':
								inshell_proteins_single_8[3] += 1
							elif residue2.get_resname() == 'CYS':
								inshell_proteins_single_8[4] += 1
							elif residue2.get_resname() == 'GLU':
								inshell_proteins_single_8[5] += 1
							elif residue2.get_resname() == 'GLN':
								inshell_proteins_single_8[6] += 1
							elif residue2.get_resname() == 'GLY':
								inshell_proteins_single_8[7] += 1
							elif residue2.get_resname() == 'HIS':
								inshell_proteins_single_8[8] += 1
							elif residue2.get_resname() == 'ILE':
								inshell_proteins_single_8[9] += 1
							elif residue2.get_resname() == 'LEU':
								inshell_proteins_single_8[10] += 1
							elif residue2.get_resname() == 'LYS':
								inshell_proteins_single_8[11] += 1
							elif residue2.get_resname() == 'MET':
								inshell_proteins_single_8[12] += 1
							elif residue2.get_resname() == 'PHE':
								inshell_proteins_single_8[13] += 1
							elif residue2.get_resname() == 'PRO':
								inshell_proteins_single_8[14] += 1
							elif residue2.get_resname() == 'SER':
								inshell_proteins_single_8[15] += 1
							elif residue2.get_resname() == 'THR':
								inshell_proteins_single_8[16] += 1
							elif residue2.get_resname() == 'TRP':
								inshell_proteins_single_8[17] += 1
							elif residue2.get_resname() == 'TYR':
								inshell_proteins_single_8[18] += 1
							elif residue2.get_resname() == 'VAL':
								inshell_proteins_single_8[19] += 1
								
			print(inshell_proteins_single_5)
			print(inshell_proteins_single_6)
			print(inshell_proteins_single_7)
			print(inshell_proteins_single_8)
		except:
			print("Error")
	except:
		print("Error2")
		pass
	inshell_proteins_5.append(inshell_proteins_single_5)
	inshell_proteins_6.append(inshell_proteins_single_6)
	inshell_proteins_7.append(inshell_proteins_single_7)
	inshell_proteins_8.append(inshell_proteins_single_8)

filename = 'NF2_5.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_5 ,outfile)
outfile.close()

filename = 'NF2_6.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_6 ,outfile)
outfile.close()

filename = 'NF2_7.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_7 ,outfile)
outfile.close()

filename = 'NF2_8.pickle'
outfile = open(filename,'wb')
pickle.dump(inshell_proteins_8 ,outfile)
outfile.close()

#EC

# Protein Class
import pandas as pd
import pickle
from sklearn.preprocessing import LabelEncoder

# dataset import and processing
ds = pd.read_csv('newdataset.csv')
pdb = list(ds.iloc[:,1].values)
new_pdb = []
for i in pdb:
	i = i.replace('.pdb', '')
	new_pdb.append(i)

pdb = new_pdb
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

func = []
for i in range(len(pdb)):
	try:
		path_ = '/content/pdb/' + str(pdb[i]).lower() + '.pdb'
		f = open(path_, "r")
	except:
		path_ = '/content/pdb/' + str(pdb[i]).upper() + '.pdb'
		f = open(path_, "r")

	for x in f:
		x = x.replace("HEADER    ", "")
		x = x.split(' ')
		ind_func = ''
		for j in range(5):
			ind_func += x[j]
		print(ind_func)
		func.append(ind_func)
		break

filename = 'NF3.pickle'
outfile = open(filename,'wb')
pickle.dump(func ,outfile)
outfile.close()

le = LabelEncoder()
func = le.fit_transform(func)
print(le.classes_, len(le.classes_))

filename = 'NF3_le.pickle'
outfile = open(filename,'wb')
pickle.dump(func ,outfile)
outfile.close()

#Motifs

# Motifs
# Takes Sequence from the PDB Structure now.
# Libraries
import pandas as pd
import numpy as np
import wget
import os
import re
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import pickle
from sklearn.preprocessing import StandardScaler, LabelEncoder
from seq_extract import get_sequence

# Tasks
# Separate window sizes (3, 5, 7, 9, 11, 13)
window = 13

# read structure from file
ds = pd.read_csv('newdataset.csv')
pdb = list(ds.iloc[:,1].values)
new_pdb = []
for i in pdb:
	i = i.replace('.pdb', '')
	new_pdb.append(i)

pdb = new_pdb
res = ds.iloc[:,2]
chain = ds.iloc[:,3]

# Identify Motifs: CC, CXC, CX4C, CX3C, CX2C, CX2CX5C, CX2CX2C, CX2CX2CX3C
motif13 = []
motif11 = []
motif9 = []
motif7 = []
motif5 = []
motif3 = []

for i in range(len(pdb)):
	print(i, len(motif13))
	sing_motif13 = np.zeros(8)
	sing_motif11 = np.zeros(8)
	sing_motif9 = np.zeros(8)
	sing_motif7 = np.zeros(8)
	sing_motif5 = np.zeros(8)
	sing_motif3 = np.zeros(8)
	
	pdb_id = str(pdb[i])
	string = ''
	print(pdb_id, i, len(pdb))
	list_ind = 0
	file1 = '/content/pdb/' + pdb_id.upper() +'.pdb'
	file2 = '/content/pdb/' + pdb_id.lower() +'.pdb'
	if os.path.isfile(file1) == True:
		string13, string11, string9, string7, string5, string3 = get_sequence(file1, res[i], chain[i], window)
	else:
		string13, string11, string9, string7, string5, string3 = get_sequence(file2, res[i], chain[i], window)
	
	# Motif check: CC, CXC, CX4C, CX3C, CX2C, CX2CX5C, CX2CX2C, CX2CX2CX3C
	if len(re.findall(r"CC", string13)) > 0:
		sing_motif13[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string13)) > 0:
		sing_motif13[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string13)) > 0:
		sing_motif13[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string13)) > 0:
		sing_motif13[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string13)) > 0:
		sing_motif13[7] = 1.0

	if len(re.findall(r"CC", string11)) > 0:
		sing_motif11[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string11)) > 0:
		sing_motif11[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string11)) > 0:
		sing_motif11[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string11)) > 0:
		sing_motif11[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string11)) > 0:
		sing_motif11[7] = 1.0


	if len(re.findall(r"CC", string9)) > 0:
		sing_motif9[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string9)) > 0:
		sing_motif9[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string9)) > 0:
		sing_motif9[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string9)) > 0:
		sing_motif9[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string9)) > 0:
		sing_motif9[7] = 1.0


	if len(re.findall(r"CC", string7)) > 0:
		sing_motif7[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string7)) > 0:
		sing_motif7[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string7)) > 0:
		sing_motif7[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string7)) > 0:
		sing_motif7[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string7)) > 0:
		sing_motif7[7] = 1.0

	if len(re.findall(r"CC", string5)) > 0:
		sing_motif5[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string5)) > 0:
		sing_motif5[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string5)) > 0:
		sing_motif5[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string5)) > 0:
		sing_motif5[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string5)) > 0:
		sing_motif5[7] = 1.0


	if len(re.findall(r"CC", string3)) > 0:
		sing_motif3[0] = 1.0
	if len(re.findall(r"C[A-Z]C", string3)) > 0:
		sing_motif3[1] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[2] = 1.0
	if len(re.findall(r"C[A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[3] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C", string3)) > 0:
		sing_motif3[4] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z][A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[5] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C", string3)) > 0:
		sing_motif3[6] = 1.0
	if len(re.findall(r"C[A-Z][A-Z]C[A-Z][A-Z]C[A-Z][A-Z][A-Z]C", string3)) > 0:
		sing_motif3[7] = 1.0

	print(sing_motif13)
	print(sing_motif11)
	print(sing_motif9)
	print(sing_motif7)
	print(sing_motif5)
	print(sing_motif3)
	motif13.append(sing_motif13)
	motif11.append(sing_motif11)
	motif9.append(sing_motif9)
	motif7.append(sing_motif7)
	motif5.append(sing_motif5)
	motif3.append(sing_motif3)

filename = 'NF4_' + str(13) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif13,outfile)
outfile.close()

filename = 'NF4_' + str(11) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif11,outfile)
outfile.close()

filename = 'NF4_' + str(9) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif9,outfile)
outfile.close()

filename = 'NF4_' + str(7) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif7,outfile)
outfile.close()

filename = 'NF4_' + str(5) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif5,outfile)
outfile.close()

filename = 'NF4_' + str(3) + '.pickle'
outfile = open(filename,'wb')
pickle.dump(motif3,outfile)
outfile.close()
