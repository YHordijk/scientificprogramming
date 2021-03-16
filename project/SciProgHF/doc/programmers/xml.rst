:orphan:
 

Documentation for Fortran XML libraries
=======================================

These libraries can read an xml document into fortran or write fortran
data to xml. For more documentation on xml see
`wikipedia <http://en.wikipedia.org/wiki/XML>`_.

::

    To improve performance a special attribute has been defined for tags:
    <tag type="realdata">. When this attribute is the data in this tag will be
    treated as a huge chunck of real data. The dimensions of this block have to be
    set by the attributes size=, which indicates the number of lines and width=,
    which indicates the number of columns.


Usage
-----

The libraries can be accesed by including in the fortran code

::

    USE xml_parser

A given XML document is read at once to memory, from where any tag can
easily

be accessed.

The next thing to be done is to declare a pointer to an XML document.

::

    type(xml_tag):: tag

Then the following functions can be called:

-  xml\_tag=>xml\_open(treename)

Starts an new xml document tree, with name treename. Returns a pointer
to the xml tree

-  xml\_tag=>xml\_read(filename,treename)

Reads a file named with filename, containing an xml document into the
document tree with name treename (name is used if treename is missing).
Returns a pointer to the xmltree.

-  call xml\_write(filename,treename)

Writes the document tree treename into an xml file, with filename.

-  xml\_tag=>xml\_close(treename)

Removes the document tree treename from memory, returns the
null-pointer.

-  logical=xml\_tag\_next(xml\_tag,name,nameattribute)

Sets the pointer xml\_tag to the next tag with name and having
optionally the attribute name to nameattribute (ie. <name
name="nameattribute" ... >). If name is "" or missing then the pointer
xml\_tag is set to the next subtag of the last non empty search name (or
the entire document). The function returns true or false depending on
wether a next tag is found. So in order to search all the subtags of say
the tag 'grid' it's sufficient to use:

::

    type(xml_tag),pointer:: tag
    !
    tag=>open('xmlfile')
    do while(xml_tag_next(tag,'grid'))
      do while(xml_tag_next(tag))
        ... code ...
      end do
    end do
    tag=>xml_close('xmlfile')

-  name=xml\_name(tag) / call xml\_name(tag,tag\_name)

returns the name of the tag or it sets the name of a tag.

-  value=xml\_attribute(tag,attribute\_name) / call
   xml\_attribute(tag,attribute\_name,value)

returns the value of an tag attribute, or sets it. xml\_name returns the
name of the tag or it sets the name of a tag.

-  real\_pointer(:,:) / char\_pointer(:) => xml\_get\_data(xml\_tag)

Returns a pointer to an array containing the character of real data that
is stored in the tag pointed to by xml\_tag. So usage could be

::

    character,pointer:: char_data(:)
    real,pointer:: real_data(:,:)
    !
    ... some code ...
    char_data=>xml_get_data(tag)
    real_data=>xml_get_data(tag)
    ... more code ...

The function returns null() if no data is stored.

-  call xml\_add\_data(xml\_tag,real / char\_data)

Writes real or character data to the current tag stored in the tag
pointed to by xml\_tag. So usage could be

::

    character(len=...),pointer:: char_data(:)
    character(len=...):: string
    real:: real_data(:,:)
    !
    ... some code ...
    call xml_add_data(tag,char_data)
    call xml_add_data(tag,string)
    call xml_add_data(tag,real_data)
    ... more code ...

If real data is entered to this function a tag will automatically be
treated as a chunck of realdata. Meaning that it wil automatically get
the attributes type=*real\_data*, size and with (ie. <tag
type="realdata" size= width= >).

-  xml\_tag=>xml\_tag\_open(xml\_tag,tagname)

Adds to the current xml\_tag a new tag with name tagname and returns a
pointer to the new tag.

-  xml\_tag=>xml\_tag\_close(xml\_tag)

Returns the a pointer to the tag that contains the current tag. to the
new tag.

Example: Writing an xml document
--------------------------------

Typical use to write an xml document: As an example we write an HTML
table:

::

    use xml_parser
    type(xml_tag),pointer:: tag
    !
    tag=>xml_open('mydoc.xml')
    tag=>xml_tag_open(tag,'table')
      tag=>xml_tag_open(tag,'tr')
        tag=>xml_tag_open(tag,'td')
          ... set up a string or data ...
          call xml_attribute(tag,'valign','20')
          call xml_attribute(tag,'textcolor','#AAAAAA')
          call xml_add_data(tag,string)
        tag=>xml_tag_close(tag)
        tag=>xml_tag_open(tag,'td')
          ... similar ...
        tag=>xml_tag_close(tag)
        ... etc ...
      tag=>xml_tag_close(tag)
    tag=>xml_tag_close(tag)
    xml_write('mydoc.xml')
    xml_close('mydoc.xml')

Example: Reading an xml document
--------------------------------

Typical use to read an xml document: Let's get all the coordinates of a
bunch of atoms:

::

    use xml_parser
    type(xml_tag),pointer:: tag
    real,pointer:: coordinates(:,:)
    character(len=...):: atom_name
    !
    tag=>xml_open('mydoc.xml')
    tag=>xml_read('mydoc.xml')
    do while(xml_get_next(tag,'atom')
      atom_name=xml_attribute(tag,'name')
      do while(xml_get_next(tag))
        if (xml_name(tag)=='coordinates') then
          coordinates=>xml_get_data(tag)
          write(*,*) 'Coordinates of ',atom_name,':',coordinates
        end if
      end do
    end do
    xml_close('mydoc.xml')

xml\_settings.h
---------------

This file is included in the fortran code and contain settings that
fine-tune the performance and layout. The variables are described inside
the file.

