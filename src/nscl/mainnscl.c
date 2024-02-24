#include <stdio.h>
#include <stdlib.h>

#include "evtnscl.h"
#include "pxitree.h"

#include "TFile.h"
#include "TTree.h"

int main(int argc, char* argv[]){
FILE *fp;
fp=fopen(argv[1],"r");

TFile *evt2root = new TFile(argv[2],"RECREATE");
TTree *Pixie16 = new TTree("Pixie16","Pixie16Tree");



Pixie16->Branch("mult",&mult,"mult/I");
Pixie16->Branch("chn",chn,"chn[mult]/B");
Pixie16->Branch("crn",crn,"crn[mult]/B");
Pixie16->Branch("sln",sln,"sln[mult]/B");
Pixie16->Branch("id",id,"id[mult]/I");
Pixie16->Branch("hlen",hlen,"hlen[mult]/I");
Pixie16->Branch("elen",elen,"elen[mult]/I");
Pixie16->Branch("trlen",trlen,"trlen[mult]/I");
Pixie16->Branch("trwlen",trwlen,"trwlen[mult]/I");
Pixie16->Branch("fcode",fcode,"fcode[mult]/I");
Pixie16->Branch("pxitime",pxitime,"pxitime[mult]/L");
Pixie16->Branch("cfdtime",cfdtime,"ctime[mult]/I");
Pixie16->Branch("cfdtimef",cfdtimef,"ctimef[mult]/I");
Pixie16->Branch("energy",energy,"energy[mult]/I");
Pixie16->Branch("extra",extra,"extra[mult]/I");
//Pixie16->Branch("tr",tr,"tr[mult][4096]/S");
//Pixie16->Branch("esum",esum,"esum[mult][4]/I");
//Pixie16->Branch("qsum",qsum,"qsum[mult][8]/I");




unsigned int sub[2016];
memset(sub,0,sizeof(sub));
struct subevent subevt[1000];
evtreader(sub,subevt, fp,Pixie16);

evt2root->cd();
Pixie16->Write();
evt2root->Close();

return 0;
}
