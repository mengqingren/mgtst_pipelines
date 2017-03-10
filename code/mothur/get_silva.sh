#!/usr/bin/sh

REFS=${1}

wget -N -P ${REFS}/ http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz
tar xvzf ${REFS}/Silva.nr_v119.tgz -C ${REFS}/

## cleanup
rm ${REFS}/README*
rm ${REFS}/Silva.nr_v119.tgz


