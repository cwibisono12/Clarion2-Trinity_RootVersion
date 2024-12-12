#include "evtnscl.h"
#include "TTree.h"
#include "TFile.h"
#include "pxitree.h"

void evt2root(int len, struct subevent *subevt,TTree *Pixie16){
int i,j,k;

j=0;
for(i=0;i<len;i++){
if(subevt[i].crn !=1) continue;
chn[j]=subevt[i].chn;
sln[j]=subevt[i].sln;
crn[j]=subevt[i].crn;
id[j]=subevt[i].id;
hlen[j]=subevt[i].hlen;
elen[j]=subevt[i].elen;
trlen[j]=subevt[i].trlen;
trwlen[j]=subevt[i].trwlen;
energy[j]=subevt[i].energy;
fcode[j]=subevt[i].fcode;
pxitime[j]=subevt[i].time;
cfdtime[j]=subevt[i].ctime;
cfdtimef[j]=subevt[i].ctimef;
extra[j]=subevt[i].extra;
j=j+1;
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
}
mult = j;
	Pixie16->Fill();	
	
}
