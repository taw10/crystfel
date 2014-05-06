#!/bin/bash

top_srcdir=$1

CRYSTFEL_BASE_VERSION=$( cf=( `grep PACKAGE_VERSION config.h` ); echo ${cf[2]} | sed -n 's/"//gp' )
sed 's/\$u\$/'${CRYSTFEL_BASE_VERSION}'/g' $top_srcdir/version.h.in > version1.tmp
command -v git > /dev/null 2>&1
if [ $? -eq 0 ]; then
  if [ -d ".git" ]; then
    git log -1 --pretty=%B | grep 'This is CrystFEL' > /dev/null
    if [ $? -eq 0 ]; then
      CRYSTFEL_GIT_COMMIT=""
    else
      CRYSTFEL_GIT_COMMIT="+"`git rev-parse HEAD`
    fi
  fi
fi
sed 's/\$e\$/'${CRYSTFEL_GIT_COMMIT}'/g' version1.tmp > version2.tmp
diff version.h version2.tmp > /dev/null
if [ $? -ne 0 ]; then
  mv version2.tmp version.h
  rm version1.tmp
else
  rm version1.tmp version2.tmp
fi
