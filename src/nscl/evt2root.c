#include "evtnscl.h"
#include "TTree.h"
#include "TFile.h"
#include "pxitree.h"

void evt2root(int len, struct subevent *subevt,TTree *Pixie16){
int i,j,k;
mult=len;
for(i=0;i<len;i++){
if(subevt[i].crn !=1) continue;
chn[i]=subevt[i].chn;
sln[i]=subevt[i].sln;
crn[i]=subevt[i].crn;
id[i]=subevt[i].id;
hlen[i]=subevt[i].hlen;
elen[i]=subevt[i].elen;
trlen[i]=subevt[i].trlen;
trwlen[i]=subevt[i].trwlen;
fcode[i]=subevt[i].fcode;
pxitime[i]=subevt[i].time;
cfdtime[i]=subevt[i].ctime;
cfdtimef[i]=subevt[i].ctimef;
energy[i]=subevt[i].energy;
extra[i]=subevt[i].extra;
//	if(hlen[i]==4 && trlen[i]==0) continue;
	/*
	k=0;
	for(j=hlen[i];j<elen[i];j++){
		tr[i][j]=subevt[i].tr[j-hlen[i]+k];
		k=k+1;
		}
	*/
	/*
	if(hlen[i]==8 || hlen[i] == 16){ 
	for(j=0;j<4;j++){
		esum[i][j]=subevt[i].esum[j];
		}
	}
	*/
	/*
	if(hlen[i] == 12 || hlen[i] == 16){
	for(j=0;j<8;j++){
		qsum[i][j]=subevt[i].qsum[j];
		}
	}
	*/
	Pixie16->Fill();	
	}
}
