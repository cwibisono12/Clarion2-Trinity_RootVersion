//Standalone Program to Check if a Point is inside Cut Files

//C. Wibisono
//06/18 '23

#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TCutG.h"

int main(int argc, char * argv[]){

TFile *f=new TFile(argv[1]);
float x=atof(argv[2]);
float y=atof(argv[3]);

TCutG* cutg21 = (TCutG*)f->Get("cutg21");
TCutG* cutg22 = (TCutG*)f->Get("cutg22");
TCutG* cutg23 = (TCutG*)f->Get("cutg23");
TCutG* cutg24 = (TCutG*)f->Get("cutg24");
TCutG* cutg25 = (TCutG*)f->Get("cutg25");
TCutG* cutg26 = (TCutG*)f->Get("cutg26");
TCutG* cutg27 = (TCutG*)f->Get("cutg27");
TCutG* cutg28 = (TCutG*)f->Get("cutg28");
TCutG* cutg29 = (TCutG*)f->Get("cutg29");
TCutG* cutg210 = (TCutG*)f->Get("cutg210");

TCutG* cutg41 = (TCutG*)f->Get("cutg41");
TCutG* cutg42 = (TCutG*)f->Get("cutg42");
TCutG* cutg43 = (TCutG*)f->Get("cutg43");
TCutG* cutg44 = (TCutG*)f->Get("cutg44");
TCutG* cutg45 = (TCutG*)f->Get("cutg45");
TCutG* cutg46 = (TCutG*)f->Get("cutg46");
TCutG* cutg47 = (TCutG*)f->Get("cutg47");
TCutG* cutg48 = (TCutG*)f->Get("cutg48");
TCutG* cutg49 = (TCutG*)f->Get("cutg49");
TCutG* cutg410 = (TCutG*)f->Get("cutg410");
TCutG* cutg411 = (TCutG*)f->Get("cutg411");
TCutG* cutg412 = (TCutG*)f->Get("cutg412");
TCutG* cutg413 = (TCutG*)f->Get("cutg413");
TCutG* cutg414 = (TCutG*)f->Get("cutg414");
TCutG* cutg415 = (TCutG*)f->Get("cutg415");
TCutG* cutg416 = (TCutG*)f->Get("cutg416");

printf("x: %f,y: %f,o: %d\n",x,y,cutg21->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg22->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg23->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg24->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg25->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg26->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg27->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg28->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg29->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg210->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg41->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg42->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg43->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg44->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg45->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg46->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg47->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg48->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg49->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg410->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg411->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg412->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg413->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg414->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg415->IsInside(x,y));
printf("x: %f,y: %f,o: %d\n",x,y,cutg416->IsInside(x,y));

return 0;
}
