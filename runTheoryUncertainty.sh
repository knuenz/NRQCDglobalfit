#CS:
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,1,0,1,0,0,0,10,70,-1.5,1.5,3,4)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,1,0,1,0,3,0,10,70,-1.5,1.5,3,4)'
#
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,1,1,1,0,0,0,10,70,-2,3.5,4,3)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,1,1,1,0,3,0,10,70,-2,3.5,4,3)'

#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(2,0,0,1,0,0,0,10,70,-1.5,1.5,3,4)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(2,0,0,1,0,3,0,10,70,-1.5,1.5,3,4)'
                                     
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(2,0,1,1,0,0,0,10,70,-2,3.5,1,2)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(2,0,1,1,0,3,0,10,70,-2,3.5,1,2)'


#Ratio 3S1/1S0 NL vs NLO:
root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,2,0,1,0,3,0,3,20,1e-2,30,1,3)'

### NLO-LO systematics:
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,2,0,1,0,0,0,2.85,70,-1.5,1.5,3,4)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,2,1,1,0,0,0,2.85,70,-2,3.5,1,2)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,2,2,1,0,0,0,2.85,70,-2,3.5,1,2)'
#root -b -q 'CalcTheoryUncertaintyGraphs.cc+(0,2,3,1,0,0,0,2.85,70,-2,3.5,1,2)'

#void CalcTheoryUncertaintyGraphs(int iModelN, int iModelD, int iMeasID=0, int nRapModel=0, int iRapModel=0, int iState=0, int iFrame=0, double pTMin=0, double pTMax=0, double ResMin=0, double ResMax=0, int LegendPos=1, int LatexPos=1){ //nState=1,2,3:Upsi1S,2S,3S
