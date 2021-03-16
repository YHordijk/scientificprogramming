##  -*- Mode: Python; python-indent: 8 -*-

#
# Dirac Input
# A graphical tool for generating Dirac molecule files
#
# Christoph Jacob <jacob@few.vu.nl>, Feb 2005
#

import os, shutil 

class DiracAtom :

	def __init__ (self, type="X", x=0.0, y=0.0, z=0.0) :
		self.type = type
		self.coords = [x, y, z]

	def to_list (self) :
		return [self.type] + self.coords

	def dump (self) :
		print "Atom: ", self.type, "xyz: ", self.coords

Elements = [ '', 
	     'H' , 'He', 'Li', 'Be', 'B' ,
	     'C' , 'N' , 'O' , 'F' , 'Ne',
	     'Na', 'Mg', 'Al', 'Si', 'P' ,
	     'S' , 'Cl', 'Ar', 'K ', 'Ca',
	     'Sc', 'Ti', 'V ', 'Cr', 'Mn',
	     'Fe', 'Co', 'Ni', 'Cu', 'Zn',
	     'Ga', 'Ge', 'As', 'Se', 'Br',
	     'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
	     'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
	     'Pd', 'Ag', 'Cd', 'In', 'Sn',
	     'Sb', 'Te', 'I' , 'Xe', 'Cs',
	     'Ba', 'La', 'Ce', 'Pr', 'Nd',
	     'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
	     'Dy', 'Ho', 'Er', 'Tm', 'Yb',
	     'Lu', 'Hf', 'Ta', 'W ', 'Re',
	     'Os', 'Ir', 'Pt', 'Au', 'Hg',
	     'Tl', 'Pb', 'Bi', 'Po', 'At',
	     'Rn', 'Fr', 'Ra', 'Ac', 'Th',
	     'Pa', 'U' , 'Np', 'Pu', 'Am',
	     'Cm', 'Bk', 'Cf', 'Es', 'Fm',
	     'Md', 'No', 'Lr', 'Rf', 'Db',
	     'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
	     'Rg']


class DiracAtomType :

	def __init__ (self, type) :
		self.type = type
		try:
			self.nuc_charge = Elements.index ( type.strip() )
		except ValueError :
			self.nuc_charge = -1

		self.basis_type = 'None'

	def set_lib_basis (self, name) :
		self.basis_type = 'Library'
		self.basis_name = name

	def is_valid (self) :
		if self.nuc_charge < 0 :
			return False
		if self.basis_type == 'None' :
			return False
		
		return True

	def is_element (self) :
		return self.type in Elements

	def basis_type_str (self) :
		if self.basis_type == 'None' :
			basis_str = 'No basis set given! '
		elif self.basis_type == 'Library' :
			basis_str = 'Library basis set: ' + self.basis_name 

		return basis_str


	def basis_set_to_str (self) :

		bset_str = ""

		if self.basis_type == 'Library' :
			bset_str = "LARGE BASIS " + self.basis_name + " \n"
		else :
			print "Something went completly wrong, unknown basis set type "
		
		return bset_str

	def dump (self) :
		print "Atomtype ", self.type
		print "  nuc charge: ", self.nuc_charge
		print "  basis type: ", self.basis_type


Pointgroup_Dict = {"C1" :" 0         ", 
		   "Ci" :" 1XYZ      ", 
		   "Cs" :" 1  Z      ", 
		   "C2" :" 1XY       ", 
		   "C2h":" 2  ZXYZ   ", 
		   "C2v":" 2 Y X     ", 
		   "D2" :" 2XY  YZ   ", 
		   "D2h":" 3  Z  Y  X"}

class DiracMoleculeFile :
	
	def __init__ (self) :
		self.title1 = ""
		self.title2 = ""
		self.charge = 0
		self.auto_symmetry = True
		self.pointgroup = "C1"

		self.atomlist = []

		self.atomtypes = []

	def save (self, filename) :
		print " Dump: "
		self.dump ()

		print " Saving to file : ", filename

		if os.path.exists (filename) :
			shutil.copyfile (filename, filename + ".bak")
		f = file(filename, "w")

		f.write ("DIRAC \n")
		f.write (self.title1 + "\n")
		f.write (self.title2 + "\n")
		f.write ("C%4d%3d" % (len(self.atomtypes), self.charge)) 

		if self.auto_symmetry :
			f.write ("           \n")
		else :
			f.write (Pointgroup_Dict[self.pointgroup]+"\n")

		for at in self.atomtypes :
			
			atoms = []
			for a in self.atomlist :
				if a.type == at.type :
					atoms.append(a)

			f.write (' %9.0f%5d \n' % (at.nuc_charge, len(atoms)) )
			
			for a in atoms :
				f.write ('%-4s %14.6f %14.6f %14.6f \n' % 
					 (a.type, a.coords[0], a.coords[1], a.coords[2]) )

			f.write ( at.basis_set_to_str () )

		f.write ("\n")
		f.close ()


	def load (self, filename) :
		pass
    
	def dump (self) :
		print " Dumping Molecule File information: "
		print "   Title1: ", self.title1
		print "   Title2: ", self.title2
		print "   Charge: ", self.charge
		if self.auto_symmetry :
			print "   Automatic Symmetry detection / Linear Symmetry : "
		else :
			print "   Symmetry point group: ", self.pointgroup

		for a in self.atomlist :
			a.dump()

		for at in self.atomtypes :
			at.dump()
