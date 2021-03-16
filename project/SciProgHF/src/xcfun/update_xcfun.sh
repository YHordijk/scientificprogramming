#!/bin/bash

echo 'Updating xcfun from official repository using Git' 
echo "Note that you have to do 'svn add' on new files manually"
echo
rm -rf ../xcuptmp/
mkdir ../xcuptmp 
cd ../xcuptmp &&
git clone git@repo.ctcc.no:xcfun.git &&
rm -rf `find xcfun -name '.git*'` &&
cp -rp xcfun/* ../xcfun &&
rm -rf ../xcuptmp/ &&
echo 'XCFun updated'
