#!/bin/csh
echo
echo Start compilation at `date`
g77 -c *.f
g77 -c *.F
set LIBS=`cernlib mathlib`
echo
echo Linking with $LIBS
echo $LIBS
g77 *.o $LIBS -oups
chmod u+x ups
echo
echo Job execution starts at `date`
#
#ln -s -f $MY_DIR/dy.out fort.8
#py2 <<EoF
#100
#1
#1
#200.
#0.8, -1.
#2.4, 4.5
#0.8
#1, 5004
#EoF
##mv py.rz /afs/cern.ch/na38/data/dypp450grvlo4.rz
#mv py.ntu $MY_DIR/dy.hbook
#rm fort.8
##
#echo
#echo  -----------------------------------------------------
#echo   py2.job has FINISHED at `date`
#echo  -----------------------------------------------------
#exit
#C!======================================================================
#  Description of the values read from the input file (STIN)
#
#' Number of events to be processed: ', NEVENT
#' Type of events (reaction) to be generated: ',IEVTYP
#  1 = Drell Yan ; 2 = c-cbar
#' Nucleon nucleon combination: ', IGEN
#  1 = pp ; 2 = pn ; 3 = np ; 4 = nn
#' Lab energy: ', ELAB
#' M hat boundaries: ', CKIN1, CKIN2
#' y_lab boundaries: ', YLBMIN, YLBMAX
#' Minimum Dimuon Mass: ', XMMIN
#' Choice of PDF (0=Pythia, 1=PDFLIB): ', IPDF
