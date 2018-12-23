#include "syspara.h"

void data_out(FILE *fp2, double t, double u[])
{

	int i;

	fprintf(fp2,"%lf ",t);
	for(i=0;i<NN;i++){
		fprintf(fp2,"%10.9lf ",u[i]);
	}
	fprintf(fp2,"\n");

}

void current(FILE *fp4, FILE *fp5, FILE *fp6, FILE *fp7, FILE *fp8, FILE *fp9, FILE *fp10, double t, double u[])
{

	out_ikr(fp4,t,u);
	out_iks(fp5,t,u);
	out_ical(fp6,t,u);
	out_inaca(fp7,t,u);
	out_inak(fp8,t,u);
	out_cicr(fp9,t,u);
	out_ikach(fp10,t,u);
	//printf("t=%lf\n",t);

}

// LTCC
void out_ical(FILE *fp6, double time, double p[])
{
	fprintf(fp6,"%e %e %e %e %e\n",time,ical.junc,ical.sl,ical.ca,ical.total);
}

// Ito
/* void out_ito (double p[])
{
	fprintf(fp11,"%lf %lf\n",time,ito.ik);
} */

// Ikr 
void out_ikr (FILE *fp4, double time, double p[])
{
	fprintf(fp4,"%lf %lf\n",time,ikr.ik);
}

// Iks
void out_iks (FILE *fp5, double time, double p[])
{
	fprintf(fp5,"%lf %lf\n",time,iks.ik);
}

// Ik1
/* void out_ik1 (double p[])
{
	fprintf(fp11,"%lf %lf\n",time,ik1.ik);
} */

// Incx
void out_inaca (FILE *fp7, double time, double p[])
{
	fprintf(fp7,"%lf %lf %lf %lf\n",time,ncx.junc,ncx.sl,ncx.j);
}

// Inak
void out_inak (FILE *fp8, double time, double p[])
{
	fprintf(fp8,"%lf %lf\n",time,inak.na);
}

// CICR 
void out_cicr (FILE *fp9, double time, double p[])
{
	fprintf(fp9,"%e %e %e %e\n",time,jrel.SRCarel,jrel.SERCA,jrel.SRleak);
}

// Ik,ACh 
void out_ikach (FILE *fp10, double time, double p[])
{
	fprintf(fp10,"%lf %lf\n",time,ikach.ik);
}

