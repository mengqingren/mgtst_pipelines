#!/usr/bin/sh
REFS=${1}

wget -N -P ${REFS}/ http://www.mothur.org/w/images/2/24/Trainset10_082014.pds.tgz
tar xvzf ${REFS}/Trainset10_082014.pds.tgz -C ${REFS}/
mv ${REFS}/trainset10_082014.pds/trainset10_082014.* ${REFS}/
rm -rf ${REFS}/trainset10_082014.pds
