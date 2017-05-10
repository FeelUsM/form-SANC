#! /bin/bash
set -e

form -D PeskinNaumov=1 -D useDirac=1 -l ComptonBorn.frm
mv ComptonBorn.log ComptonBorn-Peskin.log

form -D PeskinNaumov=1 -D useDirac=0 -l ComptonBorn.frm
mv ComptonBorn.log ComptonBorn-Peskin-notDirac.log

form -D PeskinNaumov=0 -D useDirac=1 -l ComptonBorn.frm
mv ComptonBorn.log ComptonBorn-SANC.log

form -D PeskinNaumov=0 -D useDirac=0 -l ComptonBorn.frm
mv ComptonBorn.log ComptonBorn-SANC-notDirac.log

cd original

form  -D useDirac=1 -l ComptonBorn.frm
mv ComptonBorn.log ComptonBorn-withDirac.log

form  -D useDirac=0 -l ComptonBorn.frm
mv ComptonBorn.log ComptonBorn-notDirac.log

cd ..

