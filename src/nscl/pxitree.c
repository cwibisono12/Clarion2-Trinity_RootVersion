#include "TH2.h"

int mult;
char chn[416];
char crn[416];
char sln[416];
int id[416];
int hlen[416];
int elen[416];
int trlen[416];
int trwlen[416];
unsigned int energy[416];
int fcode[416];
long long int pxitime[416];
unsigned int cfdtime[416];
char cfdtimef[416];
int extra[416];
//short tr[416][4096];
//int esum[416][4];
//int qsum[416][8];

TH2 *e_cal = new TH2I("e_cal","ID vs Energy",8192,0,8192,416,0,416);
