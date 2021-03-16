#!/usr/bin/env python

#
# from a given SCF output, generates an IHFSCC input
#
# Andre S.P. Gomes <gomes@few.vu.nl>, july 2007, originally for the TRIPOD project 

import re
import sys
import copy


###############################################################################
#
# spinor class 
#
class spinor:
	"""
	Spinor class: keeps track of spinor-specific information, such
	as the irrep (boson/fermion) it 
	"""

        def __init__ (self) :
		self.boson            = "A"
		self.fermion          = "E"
                self.eigenvalue       = 0.0  # orbital energy 
                self.alfa_composition = {}   # composition of alpha in the mulliken pop analysis and their weights
                self.beta_composition = {}   # composition of alpha in the mulliken pop analysis and their weights
                self.alfa_weight      = 0.0  # fraction of alpha character for the spinor
        	self.beta_weight      = 0.0  # fraction of beta  character for the spinor
		self.occupation       = 0.0  # occupation, from zero to one
		self.degeneracy       = 2    # degeneracy
 
        def print_compositon(self):
                print "Spinor composition (from mul.pop.):"
		print "Alpha weight:",self.alfa_weight
		print "      contributors:"
#		for key in self.alfa_composition:
#			print " "*20,key,self.alfa_composition{key}
		print "Alpha weight:",self.alfa_weight
		print "      contributors:"
#		for key in self.alfa_composition:
#			print " "*20,key,self.alfa_composition{key}

	def print_eigenvalue(self, origin=0.0):
		print "  Eigenvalue: %18.8f"%(self.eigenvalue - origin),

	def print_symmetries(self):
		print "  Irrep: %6s (%6s)"%(self.boson,self.fermion),

	def print_occupation(self):
		print "  Occup: %8.4f"%(self.occupation * self.degeneracy),


	def get_eigenvalue(self):
		return self.eigenvalue

	def set_eigenvalue(self,eigenvalue):
		self.eigenvalue = eigenvalue

	def set_symmetries(self,boson,fermion):
		self.boson   = boson
		self.fermion = fermion

	def set_occupation(self,occn):
		self.occupation = occn

	def set_degeneracy(self,degn):
		self.degeneracy = degn

# aqui tem que ser numa segunda passagem, depois que 
# ja temos a energia... ai com ela le-se o output do scf
# acha a entrada onde temos ...
#
# nada pode ser feito se nao tiver mulpop feita... tem que 
# colocar isso como condicao de erro...
	def set_composition_weights(self,alpha,beta):
		print "not possible yet"
#		self.alfa_weight = alpha
#		self.beta_weight = beta 

	def set_composition_contributions(self,alpha,beta):
		print "not possible yet"

###############################################################################
#
# electronic structure class 
# 
class molecular_electronic_structure:

	def __init__ (self):
		self.spinor_list        = []
		self.active_spinors	= []
		self.active_per_irrep   = {}
		self.active_occupied    = []
		self.active_occ_per_irrep   = {}
		self.active_virtuals    = []
		self.active_vir_per_irrep   = {}
		self.fermion_irreps     = []
		self.boson_irreps       = []
		self.homo               = 0 
		self.lumo               = 0
		self.ehmin              = 0.0
		self.ehmax              = 0.0
		self.epmin              = 0.0
		self.epmax              = 0.0

	def list_per_symmetry(self, symmetry):
		print "listing eigenvalues in symmetry",symmetry 
		for spinor in self.spinors:
			if (spinor.fermion or spinor.boson) == symmetry:
				spinor.print_energy()	

	def print_spinor_list(self,list,set_energy_origin_to_homo = False):
                if set_energy_origin_to_homo :
                        energy_origin = list[self.homo].get_eigenvalue()
                else:
                        energy_origin = 0.0
		for i, spinor in enumerate(list):
                        if set_energy_origin_to_homo :
			       print "Spinor %6d"%(i - self.homo),
                        else :
			       print "Spinor %6d"%i,
			spinor.print_eigenvalue(origin=energy_origin)
			spinor.print_occupation()
			spinor.print_symmetries()
			print ""
				
	def print_active_spinors(self):
		print "there are %5d active spinors in relccsd"%len(self.active_spinors)
		print "the fermion irreps present are ",self.fermion_irreps
		for irr in self.fermion_irreps:
			print "irrep %4s contains %5d spinors"%(irr,self.active_per_irrep[irr])
			print "                   %5d occupied"%self.active_occ_per_irrep[irr]
			print "                   %5d virtuals"%self.active_vir_per_irrep[irr]

		print "the boson irreps present are ",self.boson_irreps
		self.print_spinor_list(self.active_spinors)

	def print_all_spinors(self, set_energy_origin_to_homo = False):
		self.print_spinor_list(self.spinor_list, set_energy_origin_to_homo)

	def set_homo_lumo(self):
		for i,spinor in enumerate(self.spinor_list):
			if spinor.occupation > 0.00:
				self.homo = i 
			else:
				self.lumo = i
				break

	def check_Pspace_integrity(self,nacth,nactp):
		for i,irr in enumerate(self.fermion_irreps):
			if int(nacth[2*i]) > self.active_occ_per_irrep[irr]:
				print "Inconsistent input!, more holes requested than active occupied orbitals" 
				print "requested: %5d, maximum possible: %5d"%(int(nacth[2*i]),self.active_occ_per_irrep[irr])
				sys.exit(-1)
			elif int(nactp[2*i]) > self.active_vir_per_irrep[irr]: 
				print "Inconsistent input!, more particles requested than active virtual orbitals" 
				print "requested: %5d, maximum possible: %5d"%(int(nactp[2*i]),self.active_vir_per_irrep[irr])
				sys.exit(-1)


	def set_PmSpace_thresholds(self,nacth,npmh,nactp,npmp):
		# first we do holes, backtracking from the HOMO:
		total_ih        = 0
		total_holes     = 0
		total_particles = 0
		ehmin           = 0
		ehmax           = 0
		epmax		= 0
		epmin		= 0
		symchk_hole     = {}
		symchk_part     = {}

                Ehmin_delta_to_next = 0
                ehmax_delta_to_next = 0
                epmin_delta_to_next = 0
                epmax_delta_to_next = 0
	
		for i,irr in enumerate(self.fermion_irreps):
			total_ih    = int(npmh[2*i])  + total_ih
			total_holes = int(nacth[2*i]) + total_holes
			symchk_hole[irr] = 0

#	print "homo, total holes",self.homo,total_holes,range(self.homo,self.homo-total_holes,-1)

		for j in range(self.homo,self.homo-total_holes,-1):
			ehmin = self.spinor_list[j].eigenvalue
			ehmin_delta_to_next = abs(self.spinor_list[j-1].eigenvalue - self.spinor_list[j].eigenvalue) * 0.5 
			if j > (self.homo - total_ih):
				ehmax = self.spinor_list[j].eigenvalue
				ehmax_delta_to_next = abs(self.spinor_list[j-1].eigenvalue - self.spinor_list[j].eigenvalue) * 0.5
				symchk_hole[self.spinor_list[j].fermion] = symchk_hole[self.spinor_list[j].fermion] + 1 

#	print "final j,ehmin, ehmax, j+1,e(j+1)",j,ehmin, ehmax,j+1,self.spinor_list[j+1].eigenvalue,j-1,self.spinor_list[j-1].eigenvalue

		print "\n==================================="
		print "Intermediate hamiltonian definition"
		print "===================================\n"
		print "\n emin, emax for holes     :"
		print "             exact        :   %12.8f     %12.8f"%(ehmin,ehmax)
		print "             deltas  next :   %12.8f     %12.8f"%(ehmin_delta_to_next,ehmax_delta_to_next)
		print "             between next :   %12.8f     %12.8f"%((ehmin - ehmin_delta_to_next),(ehmax - ehmax_delta_to_next))

		for i,irr in enumerate(self.fermion_irreps):
			if int(npmh[2*i]) != symchk_hole[irr]:
				print "inconsistency in requested IH setup"
				print "asked for ",int(npmh[2*i])," pm holes in irrep",irr,"but defined",symchk_hole[irr]
				sys.exit(-1)


		total_ih        = 0
		for i,irr in enumerate(self.fermion_irreps):
			total_ih         = int(npmp[2*i])  + total_ih
			total_particles  = int(nactp[2*i]) + total_particles
			symchk_part[irr] = 0

		for j in range(self.lumo,self.lumo+total_particles,1):
			epmax = self.spinor_list[j].eigenvalue
			epmax_delta_to_next = abs(self.spinor_list[j+1].eigenvalue - self.spinor_list[j].eigenvalue) * 0.5 
			if j <= (self.homo + total_ih):
				epmin = self.spinor_list[j].eigenvalue
				epmin_delta_to_next = abs(self.spinor_list[j+1].eigenvalue - self.spinor_list[j].eigenvalue) * 0.5
				symchk_part[self.spinor_list[j].fermion] = symchk_part[self.spinor_list[j].fermion] + 1 

		print "\n emin, emax for particles :"
		print "             exact        :   %12.8f     %12.8f"%(epmin,epmax)
		print "             deltas  next :   %12.8f     %12.8f"%(epmin_delta_to_next,epmax_delta_to_next)
		print "             between next :   %12.8f     %12.8f"%((epmin + epmin_delta_to_next),(epmax + epmax_delta_to_next))
		print "\n\n"
		
		for i,irr in enumerate(self.fermion_irreps):
			if int(npmp[2*i]) != symchk_part[irr]:
				print "inconsistency in requested IH setup"
				print "asked for ",int(npmp[2*i])," pm particles in irrep",irr,"but defined",symchk_part[irr]
				sys.exit(-1)

                self.ehmin              = ehmin - ehmin_delta_to_next
                self.ehmax              = ehmax - ehmax_delta_to_next 
                self.epmin              = epmin + epmin_delta_to_next 
                self.epmax              = epmax + epmax_delta_to_next 


		input = " Dirac Input:    &CCFSPC MAXIT=100, NACTH="
                for m in nacth:
                        input = input + str(m) + ","

                input = input + " NACTP="
                for m in nactp:
                        input = input + str(m) + ","

                input = input + " FSSECT=1,1,1,1,0,0, DOIH=T &END\n"
		input = input + " Dirac Input:    &CCIH EHMIN="+str(self.ehmin)+ \
                                ", EHMAX="+str(self.ehmax)+ \
                                ", EPMIN="+str(self.epmin)+ \
                                ", EPMAX="+str(self.epmax)+ \
                                " &END\n"
		print input



	def count_spinors_in_range(self,high,low):
		for irrep in self.fermion_irreps:
			self.active_per_irrep[irrep] = 0
			self.active_occ_per_irrep[irrep] = 0
			self.active_vir_per_irrep[irrep] = 0
#                       print "#debug, irrep:",irrep,"  ",self.fermion_irreps

		for spinor in self.spinor_list:
                        print spinor.fermion
                        print self.active_per_irrep[spinor.fermion]
			if spinor.eigenvalue >= float(high) and spinor.eigenvalue <= float(low):
				self.active_spinors.append(spinor)
				self.active_per_irrep[spinor.fermion] = self.active_per_irrep[spinor.fermion] + 1 
				if (spinor.occupation > 0.00 ):
					self.active_occupied.append(spinor)	
					self.active_occ_per_irrep[spinor.fermion] = self.active_occ_per_irrep[spinor.fermion] + 1
				else:
					self.active_virtuals.append(spinor)	
					self.active_vir_per_irrep[spinor.fermion] = self.active_vir_per_irrep[spinor.fermion] + 1

				
	def exchange_spinors(self,i,j):
		tmp = copy.deepcopy(self.spinor_list[i])
		self.spinor_list[i] = copy.deepcopy(self.spinor_list[j])
		self.spinor_list[j] = copy.deepcopy(tmp)
		

	def selection_sort_SpinorEnergy(self):
		for i in range(0,len(self.spinor_list)-1):
			min  = i
			for j in range(i+1,len(self.spinor_list)):
				eig1 = self.spinor_list[min].eigenvalue
				eig2 = self.spinor_list[j].eigenvalue
				if ( eig2 < eig1 ):
					min = j 
			self.exchange_spinors(i,min)
				

	def create_spinor_list(self,scf_output,expand_degeneracies):
		self.set_spinor_list(scf_output,expand_degeneracies)
		self.selection_sort_SpinorEnergy()
		self.set_homo_lumo()


	def set_spinor_list(self,scf_output,expand_degeneracies=False):
		start_regexp_definition                        = "\s+ SCF - CYCLE"
		end_regexp_definition                          = "    gap     :\s+([-+]?\d+?\.\d*) au"
		lumo_energy_regexp_definition                  = "    E\(LUMO\) :\s+([-+]?\d+?\.\d*) au"   
		homo_energy_regexp_definition                  = "  - E\(HOMO\) :\s+([-+]?\d+?\.\d*) au"
		symmetry_definition_so_regexp_definition       = "\* Fermion symmetry ([a-zA-Z][0-9]?)([gu])?"
		symmetry_definition_spinfree_regexp_definition = "\* Boson symmetry ([a-zA-Z][0-9]?)([gu])?"
		symmetry_definition_linear_regexp_definition   = "\* Block\s+([0-9]+) in ([a-zA-Z][0-9][gu ])\:  Omega \=\s+([0-9\/]+)"
		occupied_or_not_regexp_definition              = ".+ f = \s*([0-9.]+)"
	   	eigenvalues_and_degeneracy_regexp_definition   = "\s*(-?\d*\.\d+)\s+\(\s*([0-9]+)\)"

		start                        = re.compile (r''+ start_regexp_definition                      +'', re.IGNORECASE)
		end                          = re.compile (r''+ end_regexp_definition                        +'', re.IGNORECASE)
		lumo_energy                  = re.compile (r''+ lumo_energy_regexp_definition                +'', re.IGNORECASE)
		homo_energy                  = re.compile (r''+ homo_energy_regexp_definition                +'', re.IGNORECASE)
		symmetry_definition_linear   = re.compile (r''+ symmetry_definition_linear_regexp_definition +'', re.IGNORECASE)
		symmetry_definition_spinfree = re.compile (r''+ symmetry_definition_spinfree_regexp_definition +'', re.IGNORECASE)
		symmetry_definition_so       = re.compile (r''+ symmetry_definition_so_regexp_definition +''    , re.IGNORECASE)
		occupied_or_not              = re.compile (r''+ occupied_or_not_regexp_definition            +'', re.IGNORECASE)
		eigenvalues_and_degeneracy   = re.compile (r''+ eigenvalues_and_degeneracy_regexp_definition +'', re.IGNORECASE)

		reading_scf_result = False 
		spinor_occupancy   = 0.0
                fermion_symmetry   = "E1"
                boson_symmetry     = "A"

                f = file(scf_output,'r')
                lines = f.readlines()

                for i, l in enumerate(lines) :
                        if start.match (l) :
                                reading_scf_result = True
			elif end.match(l) :
                                break

			if reading_scf_result:
				if lumo_energy.match(l):
					self.lumo = lumo_energy.match(l).group(1)

  				elif homo_energy.match(l):
					self.homo = homo_energy.match(l).group(1)

				elif occupied_or_not.match(l):
					spinor_occupancy = occupied_or_not.match(l).group(1)	

				elif symmetry_definition_so.match(l):
                                           if symmetry_definition_so.match(l).group(2) is not None:
                                                      fermion_symmetry = symmetry_definition_so.match(l).group(1)+symmetry_definition_so.match(l).group(2)
                                           else:
                                                      fermion_symmetry = symmetry_definition_so.match(l).group(1)
			   		   self.fermion_irreps.append(fermion_symmetry)

				elif symmetry_definition_spinfree.match(l):
                                        if ( symmetry_definition_spinfree.match(l).group(2) == "g" ):
                                           fermion_symmetry = "E1g"
					   boson_symmetry   = symmetry_definition_spinfree.match(l).group(1)+"g"
                                        elif ( symmetry_definition_spinfree.match(l).group(2) == "u" ):
                                           fermion_symmetry = "E1u"
					   boson_symmetry   = symmetry_definition_spinfree.match(l).group(1)+"u"
                                        else:
                                           fermion_symmetry = "E1"
					   boson_symmetry   = symmetry_definition_spinfree.match(l).group(1)

					if fermion_symmetry not in self.fermion_irreps:
						self.fermion_irreps.append(fermion_symmetry)
					if boson_symmetry not in self.boson_irreps:
						self.boson_irreps.append(boson_symmetry)

				elif symmetry_definition_linear.match(l):
					fermion_symmetry = symmetry_definition_linear.match(l).group(2)
					boson_symmetry   = symmetry_definition_linear.match(l).group(3)
					if fermion_symmetry not in self.fermion_irreps:
						self.fermion_irreps.append(fermion_symmetry)
					if boson_symmetry not in self.boson_irreps:
						self.boson_irreps.append(boson_symmetry)

				elif eigenvalues_and_degeneracy.match(l):
					for spinors in eigenvalues_and_degeneracy.finditer(l) :
						sp = spinor() 
						sp.set_symmetries(boson_symmetry,fermion_symmetry)
						sp.set_eigenvalue(float(spinors.group(1)))
						sp.set_degeneracy(int(spinors.group(2)))
						sp.set_occupation(float(spinor_occupancy))
                                                if (expand_degeneracies) :
                                                        degn = int(spinors.group(2)) / 2
                                                        if (degn > 1) :
                                                                sp.set_degeneracy(2)
                                                        for _ in range(degn) :
						                self.spinor_list.append(sp)
                                                else :
						        self.spinor_list.append(sp)

                f.close()

###############################################################################
#
# option parser
# 
def parse_cmdline():

        from optparse import OptionParser

        usage   = "make_fscc_input.py [options]"+\
                  "\n\nDescription: aide to create IHFSCC calculations with DIRAC from an SCF run"

        version = "Version 0.01"
        parser  = OptionParser(usage,version=version)


	moltra_min_default = -1.0
	moltra_max_default =  4.0

        parser.add_option("--scf_output",dest="scf_output", default=None,
                     help="file containing the results of sn SCF run (default=None)")

        parser.add_option("--process_scf_only",action="store_true",dest="process_scf_only", default=False,
                     help="whether to only read in the scf orbitals and print them, or to contine to setup of fscc input (default=False)")
        parser.add_option("--expand_degeneracies",action="store_true",dest="expand_degeneracies", default=False,
                     help="whether to count and print 2*n-fold degenerate energies n times (default=False)")
        parser.add_option("--reset_energy_scale",action="store_true",dest="reset_energy_scale", default=False,
                     help="resets orbital energies by making the HOMO the origin")

        parser.add_option("--scf_input",dest="scf_input", default=None,
                     help="input scf file, will be used as template for (IH)FSCC input (default=None)")
        parser.add_option("--fscc_input",dest="fscc_input", default=None,
                     help="name of the input fscc file to be generated (default=None)")
        parser.add_option("--ihfscc_input",dest="ihfscc_input", default=None,
                     help="name of the input ihfscc file to be generated (default=None)")

        parser.add_option("--nacth",dest="nacth", default="0,0,0,0",
                     help="number of active holes in cc (default=None)")
        parser.add_option("--npmh",dest="npmh", default="0,0,0,0",
                     help="number of active holes in the Pm space (default=None)")
        parser.add_option("--nactp",dest="nactp", default="0,0,0,0",
                     help="number of active particles in cc (default=None)")
        parser.add_option("--npmp",dest="npmp", default="0,0,0,0",
                     help="number of active particles in the Pm space (default=None)")


        parser.add_option("--moltra_min",dest="moltra_min", default=moltra_min_default,
                     help="moltra minumum (default="+str(moltra_min_default)+")")
        parser.add_option("--moltra_max",dest="moltra_max", default=moltra_max_default,
                     help="moltra maximum (default="+str(moltra_max_default)+")")

        (options, args) = parser.parse_args()

# post-process the nacth, nactp, npmh, npmp so that from the usual dirac definition, 
# we get the individual occupations for the fermion irreps

        options.nacth = options.nacth.split(',')
        options.npmh  = options.npmh.split(',')
        options.nactp = options.nactp.split(',')
        options.npmp  = options.npmp.split(',')

        return options

###############################################################################
#
# dirac input
# 
class dirac_menufile:

	def __init__(self):
		self.hamiltonian = { 'iotc' : ".IOTC", 'spinfree' : ".SPINFREE", 'gaunt' : ".GAUNT"}
		self.input       = ""

	
	def print_hamiltonian_menu(self,keywords = ['iotc','gaunt']):
		input = "**HAMILTONIAN\n"
		for key in keywords:
			if key in self.hamiltonian.keys():
				input = input + self.hamiltonian[key] + "\n"

		self.input = self.input + input

	def print_moltra_menu(self,irreps,min_energy,max_energy,tolerance=0.1):
		input = "**MOLTRA\n" + \
			".SCHEME\n" +  \
			" 6\n" + \
			".ACTIVE\n"

		for irr in range(0,irreps):
			input = input + " energies   %12.4f  %24.4f  %6.4f\n"%(min_energy,max_energy,tolerance)

		self.input = self.input + input

	def print_fscc_menu(self,nacth=[],nactp=[],sector=[1,1,1,1,0,0],ih_hole=['',''],ih_part=['',''],ih=True):
		input = "\n" + \
		        " &RELCCSD TIMING=T, IPRNT=1, INTERFACE='DIRAC6', DOSORT=F, DOENER=F, DOFSPC=T, DEBUG=F &END\n" + \
		        " &CCSORT USEOE=F &END\n"+ \
		        " &CCFSPC MAXIT=100, NACTH="

		for m in nacth:
			input = input + "%3d,"%(m)

		input = input + " NACTP="

		for m in nactp:
			input = input + "%3d,"%(m)

		input = input + " FSSECT="

		for m in sector:
			input = input + "%1d,"%(m)

		if ih:
			input = input + " DOIH=T &END\n"
			if ih_hole[0] != '' and ih_hole[1] != '' and ih_part[0] != '' and ih_part[1] != '': 
				input = input + " &CCIH EHMIN=%8.4f, EHMAX=%8.4f, EPMIN=%8.4f, EPMAX=%8.4f &END\n\n"%(ih_hole[0],ih_hole[1],ih_part[0],ih_part[1])
		else:
			input = input + " DOIH=F &END\n\n"

		self.input = self.input + input

	def print_inpfile(self,nacth=[5,5,6,6],nactp=[12,12,10,10],ih_hole=[-0.8,-0.5],ih_part=[0.21,0.365]):

#	self.print_hamiltonian_menu()
#	self.print_moltra_menu(2,-1.0,4.0)
		self.print_fscc_menu(nacth,nactp,ih_hole,ih_part,ih=True)

		print self.input



###############################################################################
#
# main program
# 
def print_fscc_summary(options):
        print "\n Summary of FSCC partition request\n"
        print " number of active holes    :",options.nacth
        print " holes included in pm      :",options.npmh
        print " number of active particles:",options.nactp
        print " particles included in pm  :",options.npmp
        print "\n"
	print "\n Summary of MOLTRA definition\n"
	print " active space from %10.4f to %10.4f"%(float(options.moltra_min),float(options.moltra_max))
#	print " encompassing %10d in fermion irrep 1"%

def main():

	options = parse_cmdline()

        datafile      = options.scf_output
        expand_degn   = options.expand_degeneracies
#       inputfile     = options.fscc_input
#       templatefile  = options.fscc_template
        inputfile     = "foo"
        templatefile  = "bar"

        print "#debug: scf,template,final",datafile,templatefile,inputfile
   
        molecule = molecular_electronic_structure()
	molecule.create_spinor_list(datafile,expand_degeneracies=expand_degn)
        molecule.print_all_spinors(set_energy_origin_to_homo=options.reset_energy_scale)
        if not options.process_scf_only:
	        print_fscc_summary(options)
	        molecule.count_spinors_in_range(options.moltra_min,options.moltra_max)
	        molecule.print_active_spinors()

	        molecule.check_Pspace_integrity(options.nacth,options.nactp)
	        molecule.set_PmSpace_thresholds(options.nacth,options.npmh,options.nactp,options.npmp)

#input = dirac_menufile()
#input.print_inpfile(options.nacth,options.npmh,[molecule.ehmin,molecule.ehmax],[molecule.epmin,molecule.epmax])

main()
