#!/usr/bin/env python

# -*- coding: UTF8 -*-
#
# Dirac Input
# A graphical tool for generating Dirac molecule files
#
# Christoph Jacob <jacob@few.vu.nl>, Feb 2005
#

import pygtk
pygtk.require('2.0')
import os
import gobject, gtk

from SimpleGladeApp import SimpleGladeApp, SimpleGladeDialog
from DiracMoleculeFile import DiracAtom, DiracAtomType, DiracMoleculeFile
from DiracBasislib import relativistic_basissets, nonrelativistic_basissets

PointgroupList = ["C1", "Ci", "Cs", "C2", "C2h", "C2v", "D2", "D2h"]


class DiracInput(SimpleGladeApp):
	def __init__(self, glade_path="diracinput.glade", root="Dirac Input", domain=None):
		SimpleGladeApp.__init__(self, glade_path, root, domain)

	def new(self):

		# The list of atoms

		self.liststore_atoms = gtk.ListStore(gobject.TYPE_STRING, gobject.TYPE_DOUBLE, 
						     gobject.TYPE_DOUBLE, gobject.TYPE_DOUBLE)

		self.liststore_atoms.connect ("row-inserted", self.on_molecule_changed)
		self.liststore_atoms.connect ("row-deleted", self.on_molecule_changed)
		self.liststore_atoms.connect ("row-changed", self.on_molecule_changed)

		self.treeview_atoms.set_model(self.liststore_atoms)

		self.atomsview_namerenderer = gtk.CellRendererText()
		self.atomsview_namerenderer.set_property('editable', True)
		self.atomsview_namerenderer.connect('edited', self.on_atomname_edited)

		self.atomsview_xrenderer = gtk.CellRendererText()
		self.atomsview_xrenderer.set_property('editable', True)
		self.atomsview_xrenderer.connect('edited', self.on_atomcoord_edited, 1)

		self.atomsview_yrenderer = gtk.CellRendererText()
		self.atomsview_yrenderer.set_property('editable', True)
		self.atomsview_yrenderer.connect('edited', self.on_atomcoord_edited, 2)

		self.atomsview_zrenderer = gtk.CellRendererText()
		self.atomsview_zrenderer.set_property('editable', True)
		self.atomsview_zrenderer.connect('edited', self.on_atomcoord_edited, 3)

		self.atomsview_namecolumn = gtk.TreeViewColumn('Atomtype', 
							       self.atomsview_namerenderer, text=0)
		self.atomsview_xcolumn = gtk.TreeViewColumn('x', self.atomsview_xrenderer, text=1)
		self.atomsview_xcolumn.set_expand (True)
		self.atomsview_ycolumn = gtk.TreeViewColumn('y', self.atomsview_yrenderer, text=2)
		self.atomsview_ycolumn.set_expand (True)
		self.atomsview_zcolumn = gtk.TreeViewColumn('z', self.atomsview_zrenderer, text=3)
		self.atomsview_zcolumn.set_expand (True)

		self.treeview_atoms.append_column (self.atomsview_namecolumn)
		self.treeview_atoms.append_column (self.atomsview_xcolumn)
		self.treeview_atoms.append_column (self.atomsview_ycolumn)
		self.treeview_atoms.append_column (self.atomsview_zcolumn)

		self.selection_atoms = self.treeview_atoms.get_selection()
		self.selection_atoms.set_mode(gtk.SELECTION_SINGLE)
		self.selection_atoms.connect("changed", self.on_selection_atoms_changed)
		self.selection_atoms.emit("changed")

		# The list of atomtypes 

		self.atomtypes = {}

		self.liststore_atomtypes = gtk.ListStore(gobject.TYPE_STRING, gobject.TYPE_STRING)

		self.liststore_atomtypes.connect ("row-inserted", self.on_molecule_changed)
		self.liststore_atomtypes.connect ("row-deleted", self.on_molecule_changed)
		self.liststore_atomtypes.connect ("row-changed", self.on_molecule_changed)

		self.treeview_atomtypes.set_model (self.liststore_atomtypes)
		
		self.atomtypesview_namecolumn = gtk.TreeViewColumn ('Atomtype', 
								    gtk.CellRendererText(), text=0)
		self.atomtypesview_basiscolumn = gtk.TreeViewColumn ('Basis set', 
								     gtk.CellRendererText(), text=1)
		self.atomtypesview_basiscolumn.set_expand (True)

		self.treeview_atomtypes.append_column (self.atomtypesview_namecolumn)
		self.treeview_atomtypes.append_column (self.atomtypesview_basiscolumn)

		self.selection_atomtypes = self.treeview_atomtypes.get_selection()
		self.selection_atomtypes.set_mode (gtk.SELECTION_SINGLE)
		self.selection_atomtypes.connect("changed", self.on_selection_atomtypes_changed)
		self.selection_atomtypes.emit("changed")
		
		# other stuff
		self.molecule_changed = False
		self.new_molecule ()


	def add_atomtype (self, atomtype) :

		self.atomtypes[atomtype.type] = atomtype

		self.liststore_atomtypes.append([atomtype.type, atomtype.basis_type_str()])

	def remove_atomtype (self, name) :

		# remove from atomtypes dictionary
		del self.atomtypes[name]

		# remove from ListStore

		iter = self.liststore_atomtypes.get_iter_first()
		while iter :
			if self.liststore_atomtypes.get_value(iter, 0) == name :
				self.liststore_atomtypes.remove(iter)
				break
			iter = self.liststore_atomtypes.iter_next (iter)

	def get_molecule (self) :

		molecule = DiracMoleculeFile ()

		molecule.title1 = self.entry_title1.get_text ()
		molecule.title2 = self.entry_title2.get_text ()
		molecule.charge = int( self.spinbutton_charge.get_value () )

		molecule.auto_symmetry = self.radiobutton_symmetry_auto.get_active ()
		if not molecule.auto_symmetry :
			index = self.combobox_pointgroup.get_active ()
			molecule.pointgroup = PointgroupList [index]

		for a in self.liststore_atoms :
			molecule.atomlist.append ( DiracAtom(a[0], a[1], a[2], a[3]) )

		for at in self.atomtypes :
			molecule.atomtypes.append ( self.atomtypes[at] )

		return molecule

	def set_molecule (self, molecule) :

		self.entry_title1.set_text (molecule.title1)
		self.entry_title2.set_text (molecule.title2)
		self.spinbutton_charge.set_value (molecule.charge)

		if molecule.auto_symmetry :
			self.radiobutton_symmetry_auto.set_active (True)
		else :
			self.radiobutton_symmetry_manual.set_active (True)
		self.combobox_pointgroup.set_active (PointgroupList.index(molecule.pointgroup))

		self.symmetry_changed ()

		self.liststore_atoms.clear ()
		for a in molecule.atomlist :
			self.liststore_atoms.append (a.to_list())

		self.selection_atoms.unselect_all ()

		self.liststore_atomtypes.clear ()
		self.atomtypes = {}
		for a in molecule.atomtypes :
			self.add_atomtype (a)
		
		self.selection_atomtypes.unselect_all ()

	def new_molecule (self) :

		if self.saved_check() :
			molecule = DiracMoleculeFile()
			self.set_molecule (molecule)
			
			self.molecule_filename = ""
			self.molecule_changed = False
			self.main_widget.set_title ("Dirac Input: [New file]")
	
	def is_saveable (self) :

		for name in self.atomtypes :
			if self.atomtypes[name].nuc_charge < 0 :
				popup = gtk.MessageDialog (self.main_widget, 0, 
							   gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
							   "No nuclear charge given \n"+
							   "for unknown Element "+ name)
				popup.run()
				popup.destroy()
				return False

			if self.atomtypes[name].basis_type == 'None' :
				popup = gtk.MessageDialog (self.main_widget, 0, 
							   gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
							   "No basis set given for atomtype " + name)

				popup.run()
				popup.destroy()
				return False

		return True

	def save_molecule (self) :

		if self.is_saveable() :
			molfile = self.get_molecule ()
			try:
				molfile.save (self.molecule_filename)
			except IOError :
				popup = gtk.MessageDialog (self.main_widget, 0, 
							   gtk.MESSAGE_ERROR, gtk.BUTTONS_OK,
							   "Error while saving file")
				popup.run()
				popup.destroy()
				return False
			else :
				self.molecule_changed = False

			return True
		else :
			return False

	def save_molecule_as (self) :

		if self.is_saveable() :
			dialog = gtk.FileChooserDialog ("Save As ...", self.main_widget, 
							gtk.FILE_CHOOSER_ACTION_SAVE,
							(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
							 gtk.STOCK_SAVE, gtk.RESPONSE_OK) )
			dialog.set_default_response(gtk.RESPONSE_OK)

			response = dialog.run()
			if (response == gtk.RESPONSE_OK) :
				self.molecule_filename = dialog.get_filename().strip()
				dialog.destroy()

				if not (self.molecule_filename[-4:] == ".mol") :
					self.molecule_filename += ".mol"

				self.main_widget.set_title("Dirac Input: " 
							   + os.path.basename(self.molecule_filename))
				return self.save_molecule ()

			else :
				dialog.destroy ()
				return False
		else:
			return False
		

	def saved_check (self) :
		
		## return value;
		##  True: everything was saved or user does not want saving
		##  False: not saved, dont quit/new whatever

		if self.molecule_changed :

			popup = gtk.MessageDialog (self.main_widget, 0, gtk.MESSAGE_QUESTION, 
						   gtk.BUTTONS_YES_NO,
						   "The current changes have not been saved. Save now?")

			popup.add_button (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
			popup.set_default_response (gtk.RESPONSE_YES)
			response = popup.run()
			popup.destroy ()

			if response == gtk.RESPONSE_YES :
				if self.molecule_filename == "" :
					return self.save_molecule_as ()
				else :
					return self.save_molecule ()
			elif response == gtk.RESPONSE_NO :
				return True
			elif response == gtk.RESPONSE_CANCEL :
				return False
		else :
			return True

	def symmetry_changed (self) :
		if self.radiobutton_symmetry_auto.get_active () :
			self.combobox_pointgroup.set_sensitive (False)
		else :
			self.combobox_pointgroup.set_sensitive (True)
		if (self.combobox_pointgroup.get_active () == -1) :
			self.combobox_pointgroup.set_active (0)

	def has_atom (self, atomname) :
		
		iter = self.liststore_atoms.get_iter_first()
		while iter :
			if self.liststore_atoms.get_value (iter, 0) == atomname :
				return True
			iter = self.liststore_atoms.iter_next (iter)
		
		return False

	## stuff that is related to the "Atoms"-ListView

	def on_selection_atoms_changed (self, widget, *args):

		(model, iter) = self.selection_atoms.get_selected ()
		if (iter == None) :
			self.button_remove_atom.set_sensitive (False)
		else :
			self.button_remove_atom.set_sensitive (True)

	def on_atomname_edited (self, widget, path, new_text) :

		old_atomname = self.liststore_atoms[path][0] 

		new_atomname = new_text[:4].strip()
		new_atomname = new_atomname[0].upper() + new_atomname[1:].lower()

		self.liststore_atoms[path][0] = new_atomname

		if not self.atomtypes.has_key(new_atomname) :
			self.add_atomtype ( DiracAtomType(new_atomname) )

		if not self.has_atom (old_atomname) :
			self.remove_atomtype (old_atomname)

		self.molecule_changed = True


	def on_atomcoord_edited (self, widget, path, new_text, column) :

		try :
			coord = float(new_text)
		except ValueError :
			pass
		else :
			self.liststore_atoms[path][column] = coord

		self.molecule_changed = True

	## stuff related to the "Atom Types"-ListView

	def on_selection_atomtypes_changed (self, widget, *args):
		pass


	def on_window_delete_event(self, widget, *args):

		if self.saved_check() :
			return gtk.FALSE
		else :
			return gtk.TRUE

	def on_new_activate(self, widget, *args):
		self.new_molecule ()

	def on_open_activate(self, widget, *args):

		## FIXME: open files !!
		print "Opening Files not implemented yet"
		pass

	def on_save_activate(self, widget, *args):
		if (self.molecule_filename == "") :
			self.save_molecule_as ()
		else :
			self.save_molecule ()

	def on_save_as_activate(self, widget, *args):
		self.save_molecule_as ()

	def on_quit_activate(self, widget, *args):
		if self.saved_check () :
			self.quit ()

	def on_about1_activate(self, widget, *args):
		about_dialog = DialogAbout (parent=self.main_widget)
		about_dialog.run ()
		about_dialog.destroy ()

	def on_molecule_changed(self, widget, *args):
		self.molecule_changed = True

	def on_symmetry_toggled(self, widget, *args):
		self.molecule_changed = True
		self.symmetry_changed ()

	def on_button_add_atom_clicked(self, widget, *args):
		self.liststore_atoms.append ( DiracAtom().to_list() )

		if not self.atomtypes.has_key("X") :
			self.add_atomtype ( DiracAtomType("X") )

	def on_button_read_atoms_clicked(self, widget, *args):
		## FIXME: reading atomic coordinates from file !!
		print "Reading atomic coordinates from file not implemented yet"

	def on_button_remove_atom_clicked(self, widget, *args):
		(model, iter) = self.selection_atoms.get_selected ()

		atomname = model.get_value (iter, 0) 
		model.remove (iter)

		if not self.has_atom (atomname) :
			self.remove_atomtype (atomname)

	def on_button_modify_atomtype_clicked(self, widget, *args):

		(model, iter) = self.selection_atomtypes.get_selected ()
		if iter :

			name = model.get_value (iter, 0) 
			path = model.get_path (iter)

			type = self.atomtypes[name] 

			dialog = DialogAtomtype (parent=self.main_widget)

			dialog.set_atomtype (type)

			response = dialog.run()
			if response == gtk.RESPONSE_OK :

				type = dialog.get_atomtype ()
				self.atomtypes[name] = type

				model[path][1] = type.basis_type_str ()

				dialog.destroy()
			else :
				dialog.destroy()


class DialogAbout(SimpleGladeDialog):
	def __init__(self, glade_path="diracinput.glade", root="dialog_about", domain=None, parent=None):
		SimpleGladeDialog.__init__(self, glade_path, root, domain, parent)

	def new(self):
		pass


class DialogAtomtype(SimpleGladeDialog):
	def __init__(self, glade_path="diracinput.glade", root="dialog_atomtype", domain=None, parent=None):
		SimpleGladeDialog.__init__(self, glade_path, root, domain, parent)

	def new(self):
		self.radiobutton_library_basis.set_active (True)
		self.radiobutton_library_basis.emit ("toggled")

		self.radiobutton_user_basis.set_active (False)
		self.radiobutton_user_basis.emit ("toggled")

		self.entry_otherbas = gtk.Entry ()
		self.entry_otherbas.show ()

		self.liststore_relbas = gtk.ListStore (gobject.TYPE_STRING)
		for i in relativistic_basissets :
			self.liststore_relbas.append ([i])

		self.combobox_relbas    = gtk.ComboBox (self.liststore_relbas)
		cell = gtk.CellRendererText()
		self.combobox_relbas.pack_start(cell, True)
		self.combobox_relbas.add_attribute(cell, 'text', 0)
		self.combobox_relbas.show ()

		self.liststore_nonrelbas = gtk.ListStore (gobject.TYPE_STRING)
		for i in nonrelativistic_basissets :
			self.liststore_nonrelbas.append ([i])

		self.combobox_nonrelbas = gtk.ComboBox (self.liststore_nonrelbas)
		cell = gtk.CellRendererText()
		self.combobox_nonrelbas.pack_start(cell, True)
		self.combobox_nonrelbas.add_attribute(cell, 'text', 0)
		self.combobox_nonrelbas.set_wrap_width(4)

		self.combobox_nonrelbas.show ()

		self.vbox_libbas.pack_start (self.combobox_relbas)

		self.combobox_libbas_oldsel = 1
		self.combobox_libbas_type.set_active (1)
		
	def set_atomtype (self, atomtype) :
		self.atomtype = atomtype

		self.main_widget.set_title ("Edit Atomtype: " + atomtype.type)
		self.label_name.set_text (atomtype.type)

		if atomtype.nuc_charge < 0 :
			self.spinbutton_charge.set_value (1)
		else :
			self.spinbutton_charge.set_value (atomtype.nuc_charge)

		if atomtype.is_element () :
			self.spinbutton_charge.set_sensitive (False)

		if atomtype.basis_type == 'Library' :

			self.radiobutton_library_basis.set_active (True)

			if atomtype.basis_name in nonrelativistic_basissets :
				self.combobox_libbas_type.set_active (0)
				i = nonrelativistic_basissets.index(atomtype.basis_name)
				self.combobox_nonrelbas.set_active  (i)
			elif atomtype.basis_name in relativistic_basissets :
				self.combobox_libbas_type.set_active (1)
				i = relativistic_basissets.index(atomtype.basis_name)
				self.combobox_relbas.set_active  (i)
			else :
				self.combobox_libbas_type.set_active (2)
				self.entry_otherbas.set_text (atomtype.basis_name)


	def get_atomtype (self) :
		
		self.atomtype.nuc_charge = self.spinbutton_charge.get_value ()

		if self.radiobutton_library_basis.get_active() :
			bastype = self.combobox_libbas_type.get_active ()
			if bastype == 0 : # non-relativistic

				iter = self.combobox_nonrelbas.get_active_iter ()
				if iter == None :
					self.atomtype.basis_type = 'None'
				else :
					basname = self.liststore_nonrelbas.get_value(iter, 0)
					self.atomtype.set_lib_basis (basname)
				
			elif bastype == 1 : # relativistic

				iter = self.combobox_relbas.get_active_iter ()
				if iter == None :
					self.atomtype.basis_type = 'None'
				else :
					basname = self.liststore_relbas.get_value(iter, 0)
					self.atomtype.set_lib_basis (basname)

			elif bastype == 2 : # Other

				basname = self.entry_otherbas.get_text ().strip()
				if basname == "" :
					self.atomtype.basis_type = 'None'
				else :
					self.atomtype.set_lib_basis (basname)
			else:
				self.atomtype.basis_type = 'None'

		else :
			print "Error getting basis set from dialog"
			self.atomtype.basis_type = 'None'


		return self.atomtype

	def on_radiobutton_library_basis_toggled(self, widget, *args):

		state = self.radiobutton_library_basis.get_active () 

		self.combobox_libbas_type.set_sensitive (state)

	def on_combobox_libbas_type_changed(self, widget, *args):

		i = self.combobox_libbas_type.get_active ()

		if self.combobox_libbas_oldsel == 0 :
			self.vbox_libbas.remove (self.combobox_nonrelbas)
		elif self.combobox_libbas_oldsel == 1 :
			self.vbox_libbas.remove (self.combobox_relbas)
		elif self.combobox_libbas_oldsel == 2 :
			self.vbox_libbas.remove (self.entry_otherbas)

		self.combobox_libbas_oldsel = i

		if i == 0 :    ## non-relativistic
			self.vbox_libbas.pack_start (self.combobox_nonrelbas)
		elif i == 1 :  ## relativistic
			self.vbox_libbas.pack_start (self.combobox_relbas)
		elif i == 2 :  ## other
			self.vbox_libbas.pack_start (self.entry_otherbas)

	def on_radiobutton_user_basis_toggled(self, widget, *args):
		pass


dirac_Input = DiracInput()
dirac_Input.run()

