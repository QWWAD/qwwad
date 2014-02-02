#! /bin/sh -e

# autogen.sh - Generate and update configuration files

test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

olddir=`pwd`
cd $srcdir

AUTORECONF=`which autoreconf`
if test -z $AUTORECONF; then
        echo "*** No autoreconf found, please install it ***"
        exit 1
fi

autoreconf -I m4 -I config --force --install

cd $olddir
test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"
