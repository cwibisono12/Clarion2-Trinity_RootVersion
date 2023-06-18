//Usage: GAGG Data Processing:
//C. Wibisono
#include "global.h"
#include "mapping.h"
#include <math.h>
//#include <stdbool.h>
#include "TCutG.h"
#include "gaggcut.h"

int gaggproc(struct sidetector *si, int parttype){
//printf("%d\n",parttype);
int gaggvalid=0;
int sicount;
int i;
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
if (parttype==1){
si[i].energy=100.*(gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}
if (parttype==2){
si[i].energy=si[i].energy=100.*(gaggquad[i-1]*pow(si[i].traceint,2.0)+gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}
}

//SiPM angle assignment:
for (sicount=0;sicount<3;sicount++){
si[i].theta[sicount]=si[i].sitheta[0][sicount];
si[i].phi[sicount]=si[i].siphi[0][sicount];
//if (sicount == 0)
//printf("theta: %f phi: %f\n",si[i].theta[sicount],si[i].phi[sicount]);
}


//Selecting GAGG events:
int ingate = 0;
    if (i == 1) { ingate = cutg21->IsInside(si[1].tpratio,si[1].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 2) { ingate = cutg22->IsInside(si[2].tpratio,si[2].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 3) { ingate = cutg23->IsInside(si[3].tpratio,si[3].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 4) { ingate = cutg24->IsInside(si[4].tpratio,si[4].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 5) { ingate = cutg25->IsInside(si[5].tpratio,si[5].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 6) { ingate = cutg26->IsInside(si[6].tpratio,si[6].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 7) { ingate = cutg27->IsInside(si[7].tpratio,si[7].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 8) { ingate = cutg28->IsInside(si[8].tpratio,si[8].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 9) { ingate = cutg29->IsInside(si[9].tpratio,si[9].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 10) { ingate = cutg210->IsInside(si[10].tpratio,si[10].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 11) { ingate = cutg41->IsInside(si[11].tpratio,si[11].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 12) { ingate = cutg42->IsInside(si[12].tpratio,si[12].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 13) { ingate = cutg43->IsInside(si[13].tpratio,si[13].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 14) { ingate = cutg44->IsInside(si[14].tpratio,si[14].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 15) { ingate = cutg45->IsInside(si[15].tpratio,si[15].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 16) { ingate = cutg46->IsInside(si[16].tpratio,si[16].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 17) { ingate = cutg47->IsInside(si[17].tpratio,si[17].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 18) { ingate = cutg48->IsInside(si[18].tpratio,si[18].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 19) { ingate = cutg49->IsInside(si[19].tpratio,si[19].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 20) { ingate = cutg410->IsInside(si[20].tpratio,si[20].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 21) { ingate = cutg411->IsInside(si[21].tpratio,si[21].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 22) { ingate = cutg412->IsInside(si[22].tpratio,si[22].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 23) { ingate = cutg413->IsInside(si[23].tpratio,si[23].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 24) { ingate = cutg414->IsInside(si[24].tpratio,si[24].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 25) { ingate = cutg415->IsInside(si[25].tpratio,si[25].traceint);
                  if (ingate ==1)gaggvalid++;}
    if (i == 26) { ingate = cutg416->IsInside(si[26].tpratio,si[26].traceint);
                  if (ingate ==1)gaggvalid++;}

//Transforming variable to another :
//Only select si event that matches condition above:


if (gaggvalid > 0 && ingate == 1){
  si[gaggvalid-1]=si[i];
}


//Time Correlation Between Germanium and GAGG:
/*
int gage;
int gev;
for (gage=1;gage<MAX_SI;gage++){
        if (si[gage].valid != 1) continue;
        si[gage].time=(si[gage].sit[0]+si[gage].sit[1])/2;
        for (gev=0;gev<gmult;gev++){
            int tdiffgage=(ge[gev].time-si[gage].time) + 2000;
            if ((ge[gev].energy > 0 && ge[gev].energy < 4096) && (tdiffgage > 0  && tdiffgage <4096)){
            gagg_ge_tdiff[ge[gev].energy][tdiffgage]++;}
                                   }     
            }

*/
 
//}//end loop over gevalid //C.W
		}// end Si condition for GAGG firing:
    
	}//end Si
return gaggvalid;
}

