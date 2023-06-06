//Standalone Program to Plot Some Spectra from ClarionTree
//C. Wibisono
//06/05 '23

#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"

#define MAX_GE 16

int main(int argc, char * argv[]){

TFile *f=new TFile(argv[1]);
TTree *T=(TTree*)f->Get("Clarion");
TFile *g=new TFile(argv[2],"RECREATE");

TH1I *forw = new TH1I("forw","",5000,0,5000);
TH1I *back = new TH1I("back","",5000,0,5000);
TH1I *perp = new TH1I("perp","",5000,0,5000);

int i,j;
int geid[MAX_GE];
int geen[MAX_GE];
char gmult;
T->SetBranchAddress("geid",geid);
T->SetBranchAddress("geenergy",geen);
T->SetBranchAddress("gemult",&gmult);

//Read Each Rows:
for(i=0;i<T->GetEntries();i++){
T->GetEvent(i);
	for(j=0;j<gmult;j++){
	if(geid[j] == 5 || geid[j] == 6 || geid[j] == 11){
	forw->Fill(geen[j]);
		}
	if(geid[j] == 1 || geid[j] == 2 || geid[j] == 7 || geid[j] == 8 || geid[j]){
	back->Fill(geen[j]); 
	}
	if(geid[j] == 3 || geid[j] == 4 || geid[j] ==9 || geid[j] ==10){
	perp->Fill(geen[j]);
	}
}
}
forw->Write();
back->Write();
perp->Write();
f->Close();
g->Close();




return 0;
}
