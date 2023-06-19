//Usage: GAGG Data Processing: --Generate PID
//C. Wibisono
#include "global.h"
#include "mapping.h"
#include <math.h>
#include <stdbool.h>
#include "gaggpid.h"
#include "TH2.h"

void gaggpid(struct gdetector *ge, struct sidetector *si, int parttype,int gmult,int Egammin, int Egammax){
int gaggvalid=0;
int sicount;
int i,k;
int gemvalid=0;
//Checking whether an event is within a gamma gate window:
for(k=0;k<gmult;k++){
if(ge[k].energy >= Egammin && ge[k].energy <= Egammax)
gemvalid++;
}
//======================================================

for (i=1;i<MAX_SI;i++){
	
  if (si[i].simult == 2 && (((si[i].siid[0] == 1 && si[i].siid[1] == 2)) || (si[i].siid[0] == 2 && si[i].siid[1] == 1)) && si[i].sipileup[0] == 0 && si[i].sipileup[1] == 0
&& (si[i].sit[0]-si[i].sit[1] + 2000 > 1990 && si[i].sit[0]-si[i].sit[1]  + 2000 < 2010 )
) 
          {
si[i].peak =  (((si[i].siqdc[0][3]-(20./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1])) + (si[i].siqdc[1][3] - (20./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/20.);
si[i].tail = (((si[i].siqdc[0][5]-(55./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1])) + (si[i].siqdc[1][5] - (55./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/20.);

si[i].traceint=(((si[i].siqdc[0][2]+si[i].siqdc[0][3]+si[i].siqdc[0][4]+si[i].siqdc[0][5]+si[i].siqdc[0][6]-(115./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1]))+(si[i].siqdc[1][2]+si[i].siqdc[1][3]+si[i].siqdc[1][4]+si[i].siqdc[1][5]+si[i].siqdc[1][6]-(115./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/30.);

si[i].tpratio=4000.*si[i].tail/si[i].peak;

si[i].id = i;
si[i].valid = 1;

if (si[i].id == gaid[i-1]){
if (parttype ==1) si[i].energy=100.*(gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
if (parttype ==2) si[i].energy=100.*(gaggquad[i-1]*pow(si[i].traceint,2.0)+gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}


//SiPM angle assignment:
for (sicount=0;sicount<3;sicount++){
si[i].theta[sicount]=si[i].sitheta[0][sicount];
si[i].phi[sicount]=si[i].siphi[0][sicount];
//if (sicount == 0)
//printf("theta: %f phi: %f\n",si[i].theta[sicount],si[i].phi[sicount]);
}



//Generating PID:
if (
    ((int)si[i].traceint > 0 && (int)si[i].traceint < 4096) && ((int)si[i].tpratio > 0 && (int)si[i].tpratio < 4096)
    && gemvalid > 0){
	
    //printf("i: %d,tpratio: %d, traceint: %d, partype: %d\n",i,(int)si[i].tpratio,(int)si[i].traceint,parttype);
	//printf("id: %d,tpratio: %d, traceint: %d\n",i,(int)si[i].tpratio,(int)si[i].traceint);
    if (i == 1){pid_qdc21->Fill((int)si[i].tpratio,(int)si[i].traceint);}
    if (i == 2){pid_qdc22->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 3){pid_qdc23->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 4){pid_qdc24->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 5){pid_qdc25->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 6){pid_qdc26->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 7){pid_qdc27->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 8){pid_qdc28->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 9){pid_qdc29->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 10){pid_qdc210->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 11){pid_qdc41->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 12){pid_qdc42->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 13){pid_qdc43->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 14){pid_qdc44->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 15){pid_qdc45->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 16){pid_qdc46->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 17){pid_qdc47->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 18){pid_qdc48->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 19){pid_qdc49->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 20){pid_qdc410->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 21){pid_qdc411->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 22){pid_qdc412->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 23){pid_qdc413->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 24){pid_qdc414->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 25){pid_qdc415->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
    if (i == 26){pid_qdc416->Fill((int)si[i].tpratio,(int)si[i].traceint);} 
}

		}// end Si condition for GAGG firing:
    
	}//end Si

}

