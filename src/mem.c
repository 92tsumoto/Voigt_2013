
#include "syspara.h"

typedef double Number;

void initial_mem()
{

	// ina_fast
	ina.Tmss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttaum=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauh=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tjss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauj=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmss==NULL || ina.Ttaum==NULL || ina.Thss==NULL || ina.Ttauh==NULL || ina.Tjss==NULL || ina.Ttauj==NULL) exit(1);
	
	// ina_late
	ina.Tmlss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauml=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thlss=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmlss==NULL || ina.Ttauml==NULL || ina.Thlss==NULL) exit(1);

	// ito
	ito.Trss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaur=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Tsss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaus=(Number *)calloc(VNMAX,sizeof(Number));
	if( ito.Trss==NULL || ito.Ttaur==NULL || ito.Tsss==NULL || ito.Ttaus==NULL) exit(1);
	
	// ikr
	ikr.Txrss=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Trkr=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikr.Txrss==NULL || ikr.Ttauxr==NULL || ikr.Trkr == NULL) exit(1);

	// iks
	iks.Txsss=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs=(Number *)calloc(VNMAX,sizeof(Number));
	if( iks.Txsss==NULL || iks.Ttauxs==NULL ) exit(1);

	// ikp
	ikp.Tss=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikp.Tss==NULL) exit(1);

	// ikur
	ikur.Txkurss=(Number *)calloc(VNMAX,sizeof(Number));
	ikur.Ttauxkur=(Number *)calloc(VNMAX,sizeof(Number));
	ikur.Tykurss=(Number *)calloc(VNMAX,sizeof(Number));
	ikur.Ttauykur=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikur.Txkurss==NULL || ikur.Ttauxkur==NULL || ikur.Tykurss == NULL || ikur.Ttauykur == NULL ) exit(1);

	// ical
	ical.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttmp1=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttmp2=(Number *)calloc(VNMAX,sizeof(Number));
	if( ical.Tdss==NULL || ical.Ttaud==NULL || ical.Tfss==NULL || ical.Ttauf==NULL || ical.Ttmp1==NULL || ical.Ttmp2==NULL) exit(1);

	// inak
	inak.Tknai=(Number *)calloc(VNMAX,sizeof(Number));
	inak.Tknao=(Number *)calloc(VNMAX,sizeof(Number));
	if( inak.Tknai==NULL || inak.Tknao==NULL ) exit(1);
	
	// inaca
	ncx.Thca=(Number *)calloc(VNMAX,sizeof(Number));
	ncx.Thna=(Number *)calloc(VNMAX,sizeof(Number));
	if( ncx.Thca==NULL || ncx.Thna==NULL ) exit(1);
	
}


void closed_mem()
{

	free(ina.Tmss); free(ina.Ttaum); free(ina.Thss); free(ina.Ttauh); free(ina.Tjss); free(ina.Ttauj);
	free(ina.Tmlss); free(ina.Ttauml); free(ina.Thlss); 
	free(ito.Trss); free(ito.Ttaur); free(ito.Tsss); free(ito.Ttaus);
	free(ikr.Txrss); free(ikr.Ttauxr); free(ikr.Trkr);
	free(iks.Txsss); free(iks.Ttauxs);
	free(ikp.Tss); 
	free(ikur.Txkurss); free(ikur.Ttauxkur); free(ikur.Tykurss); free(ikur.Ttauykur);
	free(ical.Tdss); free(ical.Ttaud); free(ical.Tfss); free(ical.Ttauf); free(ical.Ttmp1); free(ical.Ttmp2);
	free(inak.Tknai); free(inak.Tknao);
	free(ncx.Thna); free(ncx.Thca);

}

