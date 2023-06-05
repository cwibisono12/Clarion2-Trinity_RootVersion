#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "TTree.h"
#include "TFile.h"
#include "clariontree.h"

void to2root(int gaggvalid, int gmult,struct sidetector *si, struct gdetector *ge, TTree *Clarion){

	int i,j,k,l,m;
	if(gmult !=0){
	gemult=gmult;
	gaggmult=gaggvalid;
	for(k=0;k<gmult;k++){
		geid[k]=ge[k].id;
		geenergy[k]=ge[k].energy;
		geedopp[k]=ge[k].etrue;
		getime[k]=ge[k].time;
				}
	for(i=0;i<gaggvalid;i++){
		gaggid[i]=si[i].id;
		peak[i]=si[i].peak;
		tail[i]=si[i].tail;
		traceint[i]=si[i].traceint;
		gaggtime[i]=si[i].time;
		}
		Clarion->Fill();	
			}//end  loop for gmult condition


}
