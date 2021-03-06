#! /usr/bin/python

"""
Defines the classes:
- diracparser()
    Handles dirac outputfiles.
- diracnames()
   Handles the way dirac names its in and output files.


CLASS for communicating with the dirac output files.
Needs a datafile. By default the datafile: 'diracparser.dat' is used.
"""

import sys, os, os.path, time, logging

class diracparser(object):
  _keyname=[('key',''),('start',''), ('startblock', ''), ('stop', ''), ('stopblock', 0), ('type',''), ('columns',[]), ('factor',1)] 
  _log=logging.getLogger('diracparser.py')
  
  def __init__(self):
    """
    Initialization.

    _keys:     Accessible as .keys.
               Contains a list with for every key that should be read from the files.
    _keydict:  Contains the corresponding dictonaries which tells us how to read the corresponding information in the files.
    _filelist: Contains all the files we have read. Accessible with .files or with the funtions described later on.
    _energy:   Contains the eneries corresponding to these files.
               Can be read using .energy
    _datafile: Is the datafile containing information about how to read the keys.
    _resultfile: Should errormessages be written to STDOUT or to a (alread open file).
    """

    self._keys=[]
    self._keydict=[]
    self._energy=[]
    self._filelist=[]
    self._datafile='diracparser.dat'
    self._resultfile=None

  def readenergy(self, filelist):
    """
    Reads the energy.

    returns energies found in the dirac-output files in filelist.
    """
    
    energy=[]
    tmpenergy=[]
    for filename in filelist:
      if not(os.path.exists(filename)):
        if self._resultfile: self._resultfile.write('Output file: "'+filename+'" does not exist. Restart your calculation. \n')
        else: diracparser._log.error('Output file: "'+filename+'" does not exist. Restart your calculation. ')
        sys.exit(1)
      else:
        tmpdat=[]
        for key in self._keydict:
          infile=open(filename)
          lenstart=len(key['start'])
          lenblock=len(key['startblock'])
          if lenblock:
            tmp=''
            readlist=[]
            readout=0
            startblock=0
            startcol=[]
            lcount=0
            count=0
            for tmpc in infile:
              # First check for any errors
              if tmpc.count('ERROR'):
                tmpdat=['ERROR']
                break
              elif tmp>'':
                # We are in the area that we should read.
                if tmpc.count(key['stop']):
                  tmp=''
                else:
                  # We are currently in a data block.
                  if readout:
                    if (ischar and tmpc.count(key['stopblock'])) or len(tmpc)<=1:
                      readout=readout-1
                    else:
                      # Read the columns.
                      if lcount==0:
                        tmplist=[]
                        for i in range(0, len(startcol)):
                          try:
                            tmplist.append(float(key['factor'])*float(tmpc[startcol[i][0]:startcol[i][0]+startcol[i][1]]))
                          except:                            
                              if self._resultfile: self.resultfile.write('ERROR (diracparser.py): Numerical data expected. File= '+filename+'\n'+tmpc+'\n')
                              else: diracparser._log.error('Numerical data expected. File= '+filename+'\n'+tmpc)
                              tmplist.append('*')                            
                        tmplist.append('')
                        tmpenergy.append(tmplist)
                      ilast=startlen
                      count=0
                      # Get the labeling information form the block.
                      if key['type']=='label' or key['type']=='':  
                        # For type=label we read the symmetry of the orbitals. And put it at the end of the list containing the energy values.
                        for i in range(ilast, len(tmpc)-1):
                          if tmpc[i:i+1]==' ':
                            if i-ilast>1:
                              if int(tmpc[ilast:i]):
                                if len(tmpenergy[ecount][-1]): tmpenergy[ecount][-1]=tmpenergy[ecount][-1]+'<->'+readlist[lcount+count]
                                else: tmpenergy[ecount][-1]=readlist[lcount+count]
                              count=count+1
                            ilast=i
                        ecount=ecount+1
                      elif key['type'][0:5]=='label':
                        # Labeling with fixed with.
                        for i in range(ilast+1, len(tmpc), step):
                          if int(tmpc[i:i+step]):
                            if len(tmpenergy[ecount][-1]): tmpenergy[ecount][-1]=tmpenergy[ecount][-1]+'<->'+readlist[lcount+count]
                            else: tmpenergy[ecount][-1]=readlist[lcount+count]
                          count=count+1
                        ecount=ecount+1                        
                      elif key['type']=='data':  
                        # For type=data we interpret the datablock as numbers (like coordinates).
                        tmplist=[]
                        for tmpn in tmpc[ilast:].split():
                          tmplist.append(float(tmpn))
                        tmpenergy.append(tmplist)
                      else:
                        if self._resultfile: self._resultfile.write('ERROR: (diracparser.dat). Type of data block is not implemented. '+key['type']+' key= '+key['key']+' \n')
                        else: diracparser._log.error('Type of data block is not implemented. type="'+key['type']+'" key= '+key['key'])
                        sys.exit(1)
                  elif tmpc.count(key['startblock']):
                    try:
                      readout=float(key['stopblock'])
                      ischar=0
                    except:
                      readout=1
                      ischar=1
                    # If readout=0 we only need to read to the end of this line.
                    if readout==0:
                      ilast=tmpc.index(key['startblock'])+len(key['startblock'])
                      for i in range(ilast, len(tmpc)):
                        if tmpc[i:i+1] in (' ', '\n'):
                          if (i-ilast)>1:
                            try:
                              tmpdat.append(key['factor']*float(tmpc[ilast+1:i]))
                            except:                              
                              if self._resultfile: self.resultfile.write('ERROR (diracparser.py): Numerical data expected. File= '+filename+'\n'+tmpc+'\n')
                              else: diracparser._log.error('Numerical data expected. File= '+filename+'\n'+tmpc)
                              tmpdat.append('*')
                          ilast=i
                    else:
                      # Read the positions of the columns containing energy information.
                      for col in key['columns']:
                        if tmpc.count(col) and len(col)>0:
                          startcol.append((tmpc.index(col),len(col)))                      
                      startlen=tmpc.index(key['startblock'])+len(key['startblock'])
                      ilast=startlen
                      lcount+=count
                      count=0
                      ecount=0
                      # Read the labels after the block start commando.
                      if key['type'][0:5]=='label' and not(key['type']=='label'):                    
                        for i in reversed(range(0, len(key['type']))):
                          try:
                            step=int(key['type'][i:])
                          except:
                            step=int(key['type'][i+1:])
                            break
                        for i in range(ilast+1, len(tmpc)+1-step, step):
                          readlist.append(tmpc[i:i+step])
                          count=count+1
                      else:
                        for i in range(ilast, len(tmpc)-1):
                          if tmpc[i:i+1]==' ':
                            if i-ilast>1:
                              readlist.append(tmpc[ilast+1:i])
                              count=count+1
                            ilast=i
              elif tmpc.count(key['start']):
                tmp=key['start']
                tmpenergy=[]
            if tmpdat==['ERROR']:
              infile.close()
              break
            tmpdat.append(tmpenergy)                    
          else:
            # If we read only one line, find it and read all numbers on it.
            for tmpc in infile:
              if tmpc.count('ERROR'):
                tmpdat=['ERROR']
                break
              else:
                tmpn=tmpc.count(key['start'])
                if tmpn:
                  ilast=tmpc.index(key['start'])+len(key['start'])
                  if key['type'][0:5]=='label' and not(key['type']=='label'):                    
                    for i in reversed(range(0, len(key['type']))):
                      try:
                        step=int(key['type'][i:])
                      except:
                        step=int(key['type'][i+1:])
                        break
                    for i in range(ilast+1, len(tmpc)+1-step, step):
                      try:
                        tmpdat.append(key['factor']*float(tmpc[i:i+step]))
                      except:
                        if self._resultfile: self.resultfile.write('ERROR (diracparser.py): Numirical data expected. File= '+filename+'\n'+tmpc+'\n')
                        else:
                          diracparser._log.error('Numirical data expected. File= '+filename+'\n'+tmpc)                      
                  else:
                    for i in range(ilast, len(tmpc)):
                      if tmpc[i:i+1] in (' ', '\n'):
                        if (i-ilast)>1:
                          try:
                            tmpdat.append(key['factor']*float(tmpc[ilast+1:i]))
                          except:
                            if self._resultfile: self.resultfile.write('ERROR (diracparser.py): Numirical data expected. File= '+filename+'\n'+tmpc+'\n')
                            else:
                              diracparser._log.error('Numirical data expected. File= '+filename+'\n'+tmpc)
                        ilast=i
            if tmpdat==['ERROR']:
              infile.close()
              break
          infile.close()
      # Temporary debug statements
      diracparser._log.log(5,filename)
      diracparser._log.log(5,self._keys)
      diracparser._log.log(5,tmpdat)
      energy.append(tmpdat)
    return energy
  
  def set_keys(self, keys):
    """
    Makes a dictonary with how to read the output files.

    Returns nothing. For the easy to use character strings in keys the function reads the datafile.
    With this information a list of dictonaries is built and put into the variable _keydict.
    The useable names can be listed with the function keynames()
    """
    
    self._keydict=[]
    self._keys=[]

    # Connect userfriendly keywords to the characters that indicates the energy in Dirac.
    for key in keys:
      tmpdat=[]
      infile=open(self._datafile, 'r')
      if not infile:
        if self._resultfile:
          self._resultfile.write("Error (diracparser, set_keys): Datafile on how to interpret the dirac files does not exists. \n"
                       "                              File="+self._datafile)
        else:
          diracparser._log.error('Datafile on how to interpret the dirac files does not exists. \nFile='+self._datafile)
          sys.exit(1)
      for tmpc in infile:
        # Read the line starting with the current keyword into the list tmpdat.
        ilast=0
        quoted=False
        for i in range(0, len(tmpc)):
          if tmpc[i:i+1]=='#' and not(quoted): break
          elif (tmpc[i:i+1]==' ' and not(quoted)) or i==len(tmpc)-1:
            if ilast==0:
              if tmpc[ilast:ilast+1]=='"' and key==tmpc[ilast+1:i-1]:
                tmpdat=[tmpc[ilast+1:i-1]]
              elif key==tmpc[ilast:i]:
                tmpdat=[tmpc[ilast:i]]
              else:
                break
              ilast=i
            elif ilast>0:
              if tmpc[ilast:ilast+1]=='"': tmpdat.append(tmpc[ilast+1:i-1])
              else: tmpdat.append(tmpc[ilast+1:i])
            else: break
            ilast=i
          elif tmpc[i:i+1]=='"':
            quoted=not(quoted)
            if quoted: ilast=i
        if ilast: break
      else:
        if self._resultfile:
          self._resultfile.write("WARNING: Key not found in diracparser.dat. Key= ")
          self._resultfile.write(key)
        else:
          diracparser._log.error("Key not found in diracparser.dat. Key=")
          diracparser._log.error(key)
        sys.exit()
      if 0<len(tmpdat):
        # Make a dictonary from tmpdat.
        tmpkey=[]
        j=0
        for i in range(0, len(diracparser._keyname)):
          if i<len(tmpdat):
            if i==len(diracparser._keyname)-2:
              tmpkey.append((diracparser._keyname[i][0], []))
              for j in range(i, len(tmpdat)-1):
                tmpkey[i][1].append(tmpdat[j])
              j=len(tmpdat)-1
            else:          
              tmpkey.append((diracparser._keyname[i][0], tmpdat[j]))
              j=j+1
          else:
            if i==len(diracparser._keyname)-1:
              try:
                tmpn=float(tmpdat[-1])
              except:
                tmpn=0
              if (tmpn): tmpkey.append((diracparser._keyname[i][0], tmpn))
              else: tmpkey.append(diracparser._keyname[i])
            else: tmpkey.append(diracparser._keyname[i])

        self._keydict.append(dict(tmpkey))
        self._keys.append(key)
        self._energy=self.readenergy(self._filelist)

  def keynames(self):
    """
    Writes the names of the keys that could be used to STDOUT or to _resultfile.
    Returns nothing
    """
    
    infile=open(self._datafile, 'r')
    if self._resultfile: self._resultfile.write("Keys in datafile: "+self._datafile+'\n')
    else: diracparser._log.info("Keys in datafile: "+self._datafile)
    for tmpc in infile:
      quoted=False
      for i in range(0, len(tmpc)):
        if tmpc[i:i+1]=='#': break
        elif tmpc[i:i+1]==' ' and not(quoted):
          if self._resultfile: self._resultfile.write(tmpc[0:i]+'\n')
          else: print tmpc[0:i]
          break
        elif tmpc[i:i+1]=='"':
          quoted=not(quoted)
    if self._resultfile: self._resultfile.write(tmpc[0:i]+'\n')
    else: print tmpc[0:i]+'\n'
                       
  def get_keys(self):
    """
    Returns the keylist. Accessed through .keys
    """
    return self._keys

  def setfiles(self, filelist):
    """
    Set the list of files to be read Accessed through .files. Returns nothing.
    """
    self._filelist=filelist
    self._energy=self.readenergy(filelist)
    
  def addfiles(self, filelist):
    """
    Add files to the end of the filelist. Returns their read energies.
    """
    for tmpc in filelist:
      self._filelist.append(tmpc)
    tmp_energy=self.readenergy(filelist)
    for tmpdat in tmp_energy:
      self._energy.append(tmpdat)
    return tmp_energy

  def insertfile(self, pos, file):
    """
    Insert file at position pos in the file list. Returns it's energies/
    """
    self._filelist.insert(pos, file)
    tmp_energy=self.readenergy([file])
    for i in range(0, len(tmp_energy)):
      self._energy.insert(pos+i, tmp_energy[i])
    return tmp_energy[0]

  def insertfiles(self, pos, filelist):
    """
    Insert a list of files at position pos in the file list. Returns their energies.
    """
    for i in range(0, len(filelist)):
      self._filelist.insert(pos+i, filelist[i])
    tmp_energy=self.readenergy(filelist)
    for i in range(0, len(tmp_energy)):
      self._energy.insert(pos+i, tmp_energy[i])
    return tmp_energy

  def delfiles(self, filelist=[]):
    """
    Delete the named files in filelist from the _filelist.
    """
    for tmpc in filelist:
      i=index(tmpc, self._filelist)
      if i:
        del self._filelist[i]
        del self._energy[i]
      else:
        if self._resultfile:
          self._resultfile.write('WARNING (diracparser.py): File not a member of previously read files. File='+tmpc)
        else:
          diracparser._log.warning('File not a member of previously read files. File='+tmpc)

  def delfiles(self, max=0, min=0):
    """
    Delete all files between max and min from the _filelist.
    """
    # Removes the files from min up to max
    for i in range(max-1, min-1, -1):
      if i<len(self._filelist):
        del self._filelist[i]
      if i<len(self._energy):
        del self._energy[i]

  def get_files(self):
    """
    Returns the current _filelist. Accessed through .files
    """
    return self._filelist

  def get_energy(self):
    """
    Returns the current energies corresponding to the files in _filelist. Accessed through .energy
    """
    return self._energy

  def set_resultfile(self, file):
    """
    Change the resultfile. Accessed through .resultfile
    """
    self._resultfile=file

  def get_resultfile(self):
    """
    Return the current resultfile. Accessed through .resultfile
    """
    return self._resultfile

  keys=property(fget=get_keys, fset=set_keys)
  files=property(fget=get_files, fset=setfiles)
  energy=property(fget=get_energy)
  resultfile=property(fget=get_resultfile, fset=set_resultfile)

#ENDCLASS DIRACPARSER



