#!/bin/sh

tmpdir=/tmp/migrate.$$
trap "rm -fr $tmpdir" 0
trap "exit 1" 1 2 15
mkdir $tmpdir

usage()
{
  echo "usage: sh migrate.sh [--help]"                     1>&2
  echo "                     [--dir directory]"            1>&2
  echo "                     [--ignore ignored-directory]" 1>&2
  echo "                     [--ext file-extension]"       1>&2
  echo "                     [--nosave]"                   1>&2
  echo "                     [--packages-subst]"           1>&2
  echo "                     [--classes-subst]"            1>&2
  echo "                     [--dry-run]"                  1>&2
  exit $1
}

directories=
ignored=
extensions=
save=yes
packages_subst=$(dirname $0)/packages.subst
classes_subst=$(dirname $0)/classes.subst
dry_run=no
while [ $# -gt 0 ] ; do
  case "$1" in
    --help)
                 usage 0;;
    --dir)
                 shift
                 directories="$directories $1"
                 shift;;
    --ignore)
                 shift
                 ignored="$ignored $1"
                 shift;;
    --ext)
                 shift
                 extensions="$extensions $1"
                 shift;;
    --nosave)
                 save=no
                 shift;;
    --packages-subst)
                 shift
                 packages_subst="$1"
                 shift;;
    --classes-subst)
                 shift
                 classes_subst="$1"
                 shift;;
    --dry-run)
                 dry_run=yes
                 shift;;
    *)
                 echo "unknown option $1" 1>&2
                 usage 1;;
  esac
done

if [ "$directories" = "" ] ; then
  echo "missing directories" 1>&2
  usage 1
fi
if [ "$extensions" = "" ] ; then
  echo "missing extensions" 1>&2
  usage 1
fi

if [ "$packages_subst" = "" ] ; then
  echo "WARNING: package substitutions disabled at user request" 1>&2
  touch $tmpdir/migrate.sed
else
  # substitution rules for packages names
  sed -n -e 's,^ *"\([^"]*\)" *"\([^"]*\)" *,s:\1:\2:g,p' \
      < $packages_subst > $tmpdir/migrate.sed
fi
    
if [ "$classes_subst" = "" ] ; then
  echo "missing classes substitutions" 1>&2
  usage 1
fi

# substitution rules for class names, taking care of *not* substituting
# names where the pattern is only a substring of the name in the file
# so if we want to substitute OriginalClass with ReplacementClass, we
# do not want to substitute MyOwnOriginalClass or OriginalClassExtended
sed -n -e 's,^ *\([A-Za-z_][A-Za-z0-9_]*\)  *\([A-Za-z_][A-Za-z0-9_]*\) *$,s:^\1$:\2:g,p' \
    < $classes_subst >> $tmpdir/migrate.sed
sed -n -e 's,^ *\([A-Za-z_][A-Za-z0-9_]*\)  *\([A-Za-z_][A-Za-z0-9_]*\) *$,s:^\1\\([^A-Za-z0-9_]\\):\2\\1:g,p' \
    < $classes_subst >> $tmpdir/migrate.sed
sed -n -e 's,^ *\([A-Za-z_][A-Za-z0-9_]*\)  *\([A-Za-z_][A-Za-z0-9_]*\) *$,s:\\([^A-Za-z0-9_]\\)\1\$:\\1\2:g,p' \
    < $classes_subst >> $tmpdir/migrate.sed
sed -n -e 's,^ *\([A-Za-z_][A-Za-z0-9_]*\)  *\([A-Za-z_][A-Za-z0-9_]*\) *$,s:\\([^A-Za-z0-9_]\\)\1\\([^A-Za-z0-9_]\\):\\1\2\\2:g,p' \
    < $classes_subst >> $tmpdir/migrate.sed

# list the files that may be edited
find_opts=
for i in "" $ignored ; do
  if [ "$i" != "" ] ; then
    if [ "$find_opts" = "" ] ; then
      find_opts="-name $i -prune"
    else
      find_opts="$find_opts -o -name $i -prune"
    fi
  fi
done
for e in $extensions ; do
  if [ "$find_opts" = "" ] ; then
    find_opts="-name "*$e" -print"
  else
    find_opts="$find_opts -o -name "*$e" -print"
  fi
done
touch $tmpdir/files
for d in $directories ; do
  echo "find $d $find_opts"
  find $d $find_opts >> $tmpdir/files
done
if [ ! -s $tmpdir/files ] ; then
  echo "no files found!" 1>&2
  exit 1
fi

# process the files
for f in $(cat $tmpdir/files) ; do
  g=$(basename $f)
  sed -f $tmpdir/migrate.sed < $f > $tmpdir/$g
  if ! cmp -s $f $tmpdir/$g ; then
    if [ $dry_run = yes ] ; then
      echo $f
      rm $tmpdir/$g
    else
      if [ $save = yes ] ; then
        mv $f $f.orig
      fi
      mv $tmpdir/$g $f
    fi
  fi
done
