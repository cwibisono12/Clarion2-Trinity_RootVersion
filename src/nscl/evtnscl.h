//List of Declaration for Gloval Variables//
//C. Wibisono
#include "TTree.h"
#define HEADER_LENGTH 4     //unit = words with 4 bytes per word
#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

#define MAX_ID MAX_CRATES*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD


struct subevent
{
    int chn; 
    int sln;
    int crn;
    int id;
    int hlen;
    int elen;
    int trlen;          //number of samples
    int trwlen;         //number of words (two samples per word)
    int fcode;          //pileup flag
    long long int time;
    int ctime;
    int ctimef;
    int energy;
    int extra;
    short tr[4096];
    int esum[4];
    int qsum[8];        
}; 

void evtreader(unsigned int sub[], struct subevent *subevt, FILE *fpr, TTree *Pixie16);
void evt2root(int len, struct subevent *subevt, TTree *Pixie16);



