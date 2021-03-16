#!/bin/sh

#
# $Id$
#
# top-level script to generate atchives for the dyall sets
#

# preparation for the 5s/6s/7s elements
sh prepare_dyall_sblock_basis.sh

# preparation for the 4p/5p/6p elements
sh prepare_dyall_pblock_basis.sh

# preparation for the 4d/5d elements
sh prepare_dyall_dblock_basis.sh

# preparation for the 5f elements
sh prepare_dyall_fblock_basis.sh


# concatenating different blocks, first for sets with diffuses 
cat *.{4p,5p,6p,7p}.dz.aug-val.bas > dyall.av2z
cat *.{4p,5p,6p,7p}.tz.aug-val.bas > dyall.av3z
cat *.{4p,5p,6p,7p}.qz.aug-val.bas > dyall.av4z

cat *.{4p,5p,6p,7p}.dz.aug-core.bas > dyall.acv2z
cat *.{4p,5p,6p,7p}.tz.aug-core.bas > dyall.acv3z
cat *.{4p,5p,6p,7p}.qz.aug-core.bas > dyall.acv4z

# and now withouth diffuses
cat *.{4s,4p,4d,4f,5s,5p,5d,5f,6s,6p,6d,7s,7p}.dz.val.bas > dyall.v2z
cat *.{4s,4p,4d,4f,5s,5p,5d,5f,6s,6p,6d,7s,7p}.tz.val.bas > dyall.v3z
cat *.{4s,4p,4d,4f,5s,5p,5d,5f,6s,6p,6d,7s,7p}.qz.val.bas > dyall.v4z

cat *.{4s,4p,4d,4f,5s,5p,5d,5f,6s,6p,6d,7s,7p}.dz.core.bas > dyall.cv2z
cat *.{4s,4p,4d,4f,5s,5p,5d,5f,6s,6p,6d,7s,7p}.tz.core.bas > dyall.cv3z
cat *.{4s,4p,4d,4f,5s,5p,5d,5f,6s,6p,6d,7s,7p}.qz.core.bas > dyall.cv4z

# one must now edit the dyall.{a{,c}}v{2,3,4}z files and
# set the "Elements supported" section correctly by hand.
# in order to do that, one has to just collect all the 
# separate "Elements supported" fields, and leave just one
# with all the atoms
#
rm *.bas

#
# as of the time of the latest commit for this file (see Id above)
# the following can be used as template to add the list of supported
# elements to the dyall.{,c}v{2,3,4}z files above:
#
#<---- starts here
#  $
#  $ Elements Supported
#  $ K  Ca Rb Sr Cs Ba Fr Ra
#  $ Ga Ge As Se Br Kr In Sn Sb Te I  Xe Tl Pb Bi Po At Rn
#  $ Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W  Re Os Ir Pt Au Hg
#  $ Rf Db Sg Bh Hs Mt Uub Uut Uuq Uup Uuh Uus Uuo
#  $ Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr
#  $ La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
#  $
#  $  REFERENCE
#  $
#  $ 4s-7s: K.G. Dyall, J. Phys. Chem. A. (2009) 113:12638
#  $ 4p-6p: K.G. Dyall, Theor. Chem. Acc. (1998) 99:366; addendum Theor.
#  $        Chem. Acc. (2002) 108:365; revision  Theor. Chem. Acc. (2006) 115:441
#  $    4d: K.G. Dyall, Theor. Chem. Acc. (2007) 117:483.
#  $    5d: K.G. Dyall, Theor. Chem. Acc. 112, 403-409 (2004);
#  $        revision  K.G. Dyall and A.S.P. Gomes, Theor. Chem. Acc. (2009) 125:97
#  $    4f: A.S.P. Gomes, K.G. Dyall and L. Visscher DOI: 10.1007/s00214-009-0725-7
#  $    5f: K.G. Dyall, Theor. Chem. Acc. (2007) 117:491.$ 6d,7p: K.G. Dyall, in preparation
#  $ Available from the Dirac web site, http://dirac.chem.sdu.dk.
#  $
#  $
#  $
#  $
#<--- ends here
#
# to use this, remove the first 3 characters ("#  ") so that the
# dollar is at the first column (that is, use e.g in vi ":%s/#  //")
#

