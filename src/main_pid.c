//Clarion2-Trinity Software Data Processing=====================================================
//C. Wibisono
//CERN-ROOT Based Version ---06/04/2023-------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "global.h" //Global Declarations
#include "param.h" //Parameter Declarations
#include "spec_dec.h" //Macro Processor for Defining some spectra files

#include "TTree.h"
#include "TFile.h"
#include "TCutG.h"

#include "gaggcut.h"
#include "gaggpid.h"

struct subevent subevt[MAX_ID]={0};
int sevtmult=0;
unsigned long long int sevtcount=0;
unsigned long long int pileupcount=0;
unsigned long long int evtcount=0;


struct gdetector ge[MAX_GE]={0};
int gmult=0;
unsigned long long int gcount=0;


struct sidetector si[MAX_SI]={0};


int gaggvalid=0;
unsigned long long int sicount=0;
//====================================================




///////////////////////////////////
// START OF MAIN FUNCTION        //
///////////////////////////////////
int main(int argc, char **argv) {
  
  int i=0, j=0, k=0, l=0; //l is for iteration over trace
  float tempf=0;
  int max1=0, min1=0;
  int max2=0, min2=0;
  int maxid1=-1, minid1=-1;
  int maxid2=-1, minid2=-1;
  div_t e_div;
  lldiv_t lle_div;
//Parameter of kinematics Correction:
  float Ebeam,mb,mbeam,ml0,ml,mh0,mh,amu,vbz;
  int parttype=0;
  amu=931.5;
  float betascale=atof(argv[7]);
  int option=atoi(argv[8]); //use 1 to use doppler corection with kinmat or 2 to use old doppler
  int overwrite = 1;    
  int Egammin;
  int Egammax; 

  double etrace0,etrace1,btrace0,btrace1;
  double ptrace0,ptrace1,ttrace0,ttrace1,tautrace0,tautrace1;
  int dbcount = 0;
  long long int strace[500];
  memset(strace, 0, sizeof(strace));

  //temp buffer for each sub event
  unsigned int sub[MAX_SUB_LENGTH];
  memset(sub, 0, sizeof(sub));
  
  //Reference time and difference for event building
  long long int etime, tdif, idtime[MAX_ID]={0}, temptime;

  // Check that the corrent number of arguments were provided.
  if (argc<3)    {
    printf("Incorrect number of arguments:\n");
    printf("%s -op datafile calibrationfile mapfile \n", argv[0]);
    printf("\n .... calibration file is optional\n");
    printf(" .... map file is optional\n");
    printf(" .... o for overwrite spectra\n");
    printf(" .... u for update spectra\n");
    printf(" .... p for print realtime stats\n");
    return 1;
  }
    
  if(strstr(argv[1], "u") != NULL) {
    overwrite = 0;
    printf("Updating Spectra\n");
  }
  else {
    printf("Overwriting Spectra\n");
  }  
  
  //open list-mode data file from PXI digitizer  
  FILE *fpr;
  long int fprsize,fprpos;
  if ((fpr = fopen(argv[2], "r")) == NULL) {
    fprintf(stderr, "Error, cannot open input file %s\n", argv[2]);
    return 1;
  }
 
  //get file size
  fseek(fpr, 0L, SEEK_END);
  fprsize = ftell(fpr);
  rewind(fpr);
  
   
  //open debug file for streaming an 1d array
  FILE *debugfile;
  if ((debugfile = fopen(DEBUGFN, "w")) == NULL) {
    fprintf(stderr, "Error, cannot open %s\n", DEBUGFN);
    return 1;
  }   
   

//File to write GAGG_PID:
 TFile *fpidqdc = new TFile(argv[11],"RECREATE");
  char line[LINE_LENGTH];


  //open energy and time calibration file (e.g., *.ca3 file)  
  int calid=0; 
  float caloffset=0.0, calgain=0.0;
  int firstcal=0, necal=0, ntcal=0;
  
  FILE *fprcal;
  
  for (i=0; i<MAX_ID; i++) {
    ecal[i][0] = caloffset;
    ecal[i][1] = calgain;
    tcal[i][0] = caloffset;
    tcal[i][1] = calgain;     
  }  
  
  if (argc >= 4) {

    if ((fprcal = fopen(argv[3], "r")) == NULL) {
        fprintf(stderr, "Error, cannot open input file %s\n", argv[3]);
        return 1;
    } 
    
    printf("%s loaded!\n", argv[3]);

	while(fgets(line, LINE_LENGTH, fprcal) != NULL){
	    calid=0; caloffset=0; calgain=0;
        for(i=0; i<LINE_LENGTH; i++){
            if(line[i]=='#'){
                if(PRINT_CAL)printf("%s", line);
                break;
            }
            else if(line[i]>=0){
                if (firstcal==0) {
                    sscanf(line,"%f\n", &new_gain);    
	                if(PRINT_CAL) printf("%f\n", new_gain);                                       
                    firstcal=1;
                    break;
                }    
	            sscanf(line,"%d\t%f\t%f\n", &calid, &caloffset, &calgain);
	            if(PRINT_CAL) printf("%d\t%.4f\t%.4f\n", calid, caloffset, calgain);
		        if(calid >=0 && calid < MAX_ID) {
                    ecal[calid][0] = caloffset;
                    ecal[calid][1] = calgain;
                    necal++;
          		    break;
			    }
		        if(calid >=1000 && calid < 1000+MAX_ID) {
                    tcal[calid-1000][0] = caloffset;
                    tcal[calid-1000][1] = calgain;
                    ntcal++;
          		    break;
			    }			    
			    else {
				    printf("Error in reading %s : bad id or format\n", argv[3]);
				    return -1;
			    }
            }
	        else if(line[i]=='\n'){
                if(PRINT_CAL) printf("\n");
                break;
            }
            else {
                continue;
            }

        }
        	memset(line, 0, LINE_LENGTH);
 	  }
 	  
 	fclose(fprcal);
    printf("read %d energy calibrations\n", necal);
    printf("read %d time calibrations\n", ntcal);
  }
  


  //open ID->Detector map file (e.g., *.map file)  
  int mapid=0, detid=0, detidi; 
  char dettype=0;
  float theta=0, phi=0, thetai=0, phii=0, theta1=0, phi1=0, theta2=0, phi2=0, beta=0; 
  FILE *fprmap;
  int nmmap=0;
  if (argc >= 5) {

    if ((fprmap = fopen(argv[4], "r")) == NULL) {
        fprintf(stderr, "Error, cannot open input file %s\n", argv[4]);
        return 1;
    } 

    printf("%s loaded!\n", argv[4]);

	while(fgets(line, LINE_LENGTH, fprmap) != NULL){
	    mapid=0; dettype=0; thetai=0; phii=0; theta=0; phi=0; theta1=0; phi1=0; theta2=0; phi2=0; 
        for(i=0; i<LINE_LENGTH; i++){
            if(line[i]=='#'){
                if(PRINT_MAP)printf("%s", line);
                break;
            }
            else if(line[i]>=0){   
	            sscanf(line,"%d\t%c\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &mapid, &dettype, &detid, &detidi,&beta, &theta, &phi, &thetai, &phii, &theta1, &phi1, &theta2, &phi2);
	            if(PRINT_MAP) printf("%d\t%c\t%d\t%d\t%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", mapid, dettype, detid, detidi, beta, theta, phi, thetai, phii, theta1, phi1, theta2, phi2);
		        if(mapid >=0 && mapid < MAX_ID) {
                    map2type[mapid] = dettype;
                    map2det[mapid] = detid;
                    map2deti[mapid] = detidi;
		    map2beta[mapid] = beta;
                    mapangles[mapid][0] = theta;
                    mapangles[mapid][1] = phi;
                    mapanglesi[mapid][0]=thetai; 
                    mapanglesi[mapid][1]=phii; 
                    mapangles1[mapid][0]= theta1; 
                    mapangles1[mapid][1]= phi1; 
                    mapangles2[mapid][0]= theta2; 
                    mapangles2[mapid][1]= phi2; 
                    nmmap++;
          		    break;
			    }			    
			    else {
				    printf("Error in reading %s : bad id or format\n", argv[4]);
				    return -1;
			    }
            }
	        else if(line[i]=='\n'){
                if(PRINT_MAP) printf("\n");
                break;
            }
            else {
                continue;
            }

        }
        	memset(line, 0, LINE_LENGTH);
 	  }
 	  
 	fclose(fprmap);
    printf("read %d id maps\n", nmmap);
  }


//GAGG_Calibration_Parameter:


if (argc >=6){
FILE *fgagg, *fparticle, *fgamgate;
fgagg=fopen(argv[5],"r");
int idgagg;
int n;
float quadgagg,slopegagg,interceptgagg;
//int parttype; //define particle type
char type[20];
int part_type;

fparticle=fopen(argv[9],"r");
fgamgate =fopen(argv[10],"r");

if (fparticle == NULL){
fprintf(stderr,"Error, cannot open input file %s\n",argv[9]);
}

while(fgets(line,LINE_LENGTH,fparticle) !=NULL){
        if(strchr(line,'#') ==NULL)
         sscanf(line,"%f\t%f\t%f\t%f\t%s\t%d\n",&Ebeam,&mb,&ml0,&mh0,type,&parttype);
                        }

memset(line,0,LINE_LENGTH);

mbeam=mb*amu; //conversion to amu
ml=ml0*amu;
mh=mh0*amu;
vbz=pow(2.*Ebeam/mbeam,0.5);

if (fgagg == NULL){
fprintf(stderr,"Error, cannot open input file %s\n",argv[5]);
}

//Protons:
if(part_type==1){
while (fgets(line,LINE_LENGTH,fgagg) != NULL){
sscanf(line,"%d\t%f\t%f\n",&idgagg,&slopegagg,&interceptgagg);
gaid[n]=idgagg;
gaggslope[n]=slopegagg;
gaggintercept[n]=interceptgagg;
n++;
memset(line,0,LINE_LENGTH);
	}
}

//Alphas:
if(part_type==2){
while (fgets(line,LINE_LENGTH,fgagg) != NULL){
sscanf(line,"%d\t%f\t%f\t%f\n",&idgagg,&interceptgagg,&slopegagg,&quadgagg);
gaid[n]=idgagg;
gaggquad[n]=quadgagg;
gaggslope[n]=slopegagg;
gaggintercept[n]=interceptgagg;
n++;
memset(line,0,LINE_LENGTH);
		}
}

//Reading The lower and upper limit of Gamma Energy for making Banana Gate
if (fgamgate == NULL){
fprintf(stderr,"Error, cannot open input file %s\n",argv[5]);
}

while(fgets(line,LINE_LENGTH,fgamgate) !=NULL){
        if(strchr(line,'#') ==NULL)
                sscanf(line,"%d\t%d\n",&Egammin,&Egammax);
                        }
memset(line,0,LINE_LENGTH);

printf("Emin: %d,Emax: %d\n",Egammin,Egammax);
fclose(fgagg);
fclose(fparticle);
fclose(fgamgate);
}

//Open root cutg files: \\C.W
if(argc >= 7){
TFile *fcutg = new TFile(argv[6]); 
cutg21 = (TCutG*)fcutg->Get("cutg21");
cutg22 = (TCutG*)fcutg->Get("cutg22");
cutg23 = (TCutG*)fcutg->Get("cutg23");
cutg24 = (TCutG*)fcutg->Get("cutg24");
cutg25 = (TCutG*)fcutg->Get("cutg25");
cutg26 = (TCutG*)fcutg->Get("cutg26");
cutg27 = (TCutG*)fcutg->Get("cutg27");
cutg28 = (TCutG*)fcutg->Get("cutg28");
cutg29 = (TCutG*)fcutg->Get("cutg29");
cutg210 = (TCutG*)fcutg->Get("cutg210");
cutg41 = (TCutG*)fcutg->Get("cutg41");
cutg42 = (TCutG*)fcutg->Get("cutg42");
cutg43 = (TCutG*)fcutg->Get("cutg43");
cutg44 = (TCutG*)fcutg->Get("cutg44");
cutg45 = (TCutG*)fcutg->Get("cutg45");
cutg46 = (TCutG*)fcutg->Get("cutg46");
cutg47 = (TCutG*)fcutg->Get("cutg47");
cutg48 = (TCutG*)fcutg->Get("cutg48");
cutg49 = (TCutG*)fcutg->Get("cutg49");
cutg410 = (TCutG*)fcutg->Get("cutg410");
cutg411 = (TCutG*)fcutg->Get("cutg411");
cutg412 = (TCutG*)fcutg->Get("cutg412");
cutg413 = (TCutG*)fcutg->Get("cutg413");
cutg414 = (TCutG*)fcutg->Get("cutg414");
cutg415 = (TCutG*)fcutg->Get("cutg415");
cutg416 = (TCutG*)fcutg->Get("cutg416");

fcutg->Close();
}
//=============================================//
  
  /////////////////////
  // MAIN WHILE LOOP //
  /////////////////////
  while (1) { //main while loop 
  
  //Read to file and build event; //C.W
  sevtmult=toreader(sub, subevt, fpr);
     
  if (sevtmult==0) break; //end main WHILE LOOP when out of events 
  mult[0][sevtmult]++; //Histogram raw sub event multiplicity 
  sevtcount += sevtmult;
  evtcount++; //event-built number
      /////////////////////////////////////
      // END UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////////


      //skip detector building below if no map file
      if (argc >= 5) {        


      memset(&ge, 0, sizeof(ge)); //This is needed but could be replaced by setting suppress, pileup, nonprompt, clean, x/s/bmult to zero at start of loop!
      memset(&si,0,sizeof(si)); //C. W

//Start Mapping Raw Subevents into Detector Types:
detmaps(sevtmult, subevt, ge, si);


//Count Germanium multiplicities:
gmult=gamproc(ge);

//Create GAGG PID:
gaggpid(ge,si,parttype,gmult,Egammin,Egammax);

      } //end argc >= 5 condition 

      //event stats, print status every 10000 events
      lle_div=lldiv(evtcount,10000);
      if ( lle_div.rem == 0 && strstr(argv[1], "p") != NULL) {
        fprpos = ftell(fpr);
        tempf = (float)fprsize/(1024.*1024.*1024.);
        printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\r", sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
      }      

  
//printf("Ep<200: %d,Ep>200: %d\n",plt200,pgt200);
        
  } // end main while loop 
  /////////////////////////
  // END MAIN WHILE LOOP //
  /////////////////////////
  fprpos = ftell(fpr);
  tempf = (float)fprsize/(1024.*1024.*1024.);
  printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\r", sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
           

  
  
  
  
  
  printf("\n\n\n\nWriting...\n");
  fpidqdc->cd();
  pid_qdc21->Write();
  pid_qdc22->Write();
  pid_qdc23->Write();
  pid_qdc24->Write();
  pid_qdc25->Write();
  pid_qdc26->Write();
  pid_qdc27->Write();
  pid_qdc28->Write();
  pid_qdc29->Write();
  pid_qdc210->Write();
  pid_qdc41->Write();
  pid_qdc42->Write();
  pid_qdc43->Write();
  pid_qdc44->Write();
  pid_qdc45->Write();
  pid_qdc46->Write();
  pid_qdc47->Write();
  pid_qdc48->Write();
  pid_qdc49->Write();
  pid_qdc410->Write();
  pid_qdc411->Write();
  pid_qdc412->Write();
  pid_qdc413->Write();
  pid_qdc414->Write();
  pid_qdc415->Write();
  pid_qdc416->Write();
 
  fclose(fpr);
  fclose(debugfile);
  fpidqdc->Close();
  return 0;
}






