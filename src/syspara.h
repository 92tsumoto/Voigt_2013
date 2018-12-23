//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "/home/tsumoto/lib/xhplot.h"

#define NN 41
#define BUF 100
#define NUM 60

//#define R 8314.472		// J/mmol/K
//#define F 96485.33771638995	// C/mol
//#define T 310.0		// K
#define R 8314.0	// J/mol/K
#define F 96485.0	// C/mmol
#define T 310.0		// K

#define dvm 5
#define Emax 2000
#define Emin -2000
#define VNMAX (Emax-Emin)*dvm+1

struct varstruct {

    int datas;
    int line_wid[NUM];
	
	int n;
    double Istim;
	double dIstim;

	// An invariant constant
	double RTonF,RTon2F,RTon4F;

	// Cell tupe
	// 0: control, 1: isAF
	int celltype;
	// simulation type
	// 0: without, 1: with ISO stimulation, 2: with Ach stimulation
	int simtype;

	// Cell Geometry
	double length,a;
	double vcell,vmyo,vsr,vsl,vjunc;
	double ageo,acap,Cmem;
	double Fjunc,Fsl,Fjunc_CaL,Fsl_CaL;
	double J_na_juncsl,J_na_slmyo;
	double J_ca_juncsl,J_ca_slmyo;

	// Ion Valences 
	double zna,zk,zca;

	// time ratio
	double ndis;	

	// Reversal potential
	double Ena_junc,Ena_sl;
	double Ek,Eks;
	double Eca_junc,Eca_sl;
	double Ecl;
	double prnak;
			
	// Total Ion currents 
	double Ina_total;
	double Ik_total;
	double Ica_total;
	double Itotal;

	// leak and xfer
	double vleak, vxfer;

	// Extracellular ion concentrations
	double nao,ko,cao,clo;
	double mgi,cli;

	// carbachol (concentration; microM)
	double CCh;

	// Base Currnt Stimulus
	double Istim_base;

	// test variable
	double dt,dvdt;

	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

    int m;
    int l;

    double x0[NUM][NN];
    double tsign[NUM];
    double tend[NUM];

    int pflag;
    int write, graph;
    int write0;
    int half;
	int deb;
	int pswitch, sswitch;
	int out_data;

} var;

// potassium current
struct ikstruct {
	double tot;
} ik;

// calcium current
struct icastruct {
	double tot_junc,tot_sl;
} ica;

// cloride current
struct iclstruct {
	double tot;
} icl;

// Fast and Late sodium currnets
struct inastruct {

	double Gna,Gnal;
	double fast_junc,fast_sl,na;
	double *Tmss,*Ttaum,*Thss,*Ttauh,*Tjss,*Ttauj;
	double mss,taum,hss,tauh,jss,tauj;
	double late_junc,late_sl,late_na;
	double *Tmlss,*Ttauml,*Thlss;
	double mlss,tauml,hlss,tauhl;
	double tot_junc,tot_junc2,tot_sl,tot_sl2;

} ina;

// Transient Outward Current (Ito)
struct itostruct {

	double Gtof;
	double rss,taur,sss,taus;
	double *Trss,*Ttaur,*Tsss,*Ttaus;
	double fast,slow,ik;

} ito;

// Ultra Rapid activating potassium current (Ikur)
struct ikurstruct {

	double ik,Gkur;
	double xkurss,tauxkur,ykurss,tauykur;
	double *Txkurss,*Ttauxkur,*Tykurss,*Ttauykur;

} ikur;

// Rapid activating potassium current (Ikr)
struct ikrstruct {

	double ik,Gkr,gkr;
	double xr,xrss,tauxr,rkr;
	double *Txrss,*Ttauxr,*Trkr;

} ikr;

// Slowlactivating potassium current (Iks)
struct iksstruct {

	double Gks,gks_junc,gks_sl;
	double xsss,tauxs;
	double *Txsss,*Ttauxs;
	double junc,sl,ik;
		
} iks;

// plateau K Current (Ikp)
struct ikpstruct {
	double Gkp,junc,sl,ik;
	double ss,*Tss;

} ikp;

// Inward rectifier potassium current (Ik1)
struct ik1struct {

	double ik,Gk1max,gk1;
	double k1ss,ak1,bk1;

} ik1;

// Ca-activated Cl current (ICaCl)
struct iclcastruct {

	double Gclca,kd_cl_ca;
	double junc,sl,cl;

} iclca;

// L-type Calcium channel current (IcaL)
struct icalstruct {

	double dss,taud,fss,tauf;
	double *Tdss,*Ttaud,*Tfss,*Ttauf;
	double tmp1,tmp2,*Ttmp1,*Ttmp2;
	double fcaCaMSL,fcaCaj,fcaCa;
	double bar_ca_j,bar_ca_sl;
	double bar_k;
	double bar_na_j,bar_na_sl;
	double junc,sl,ca;
	double k;
	double cana_j,cana_sl,cana;
	double pca,pk,pna,Q10CaL,Qpow,total;

} ical;

// Na-K Pump
struct inakstruct {

	double kmna,kmk,Q10nak,Q10kmnai;
	double knai,knao,*Tknai,*Tknao;
	double sigma,fnak,Gmax;
	double junc,sl,na;

} inak;

// Na-Ca exchanger
struct ncxstruct {

	double hca,hna;
	double *Thca,*Thna;
	double s1_junc,s1_sl;
	double s2_junc,s2_sl;
	double s3_junc,s3_sl;
	double ksat,gamma;
	double Ka_junc,Ka_sl;
	double kmnai,kmnao,kmcai,kmcao,kdact;
	double Gmax,Q10NCX;
	double junc,sl,j;

} ncx;

// Sarcolemmal Ca Pump
struct ipcastruct {

	double Q10SLCaP;
	double Pmax,kmPca;
	double junc,sl,ca;

} ipca;

// Na Background Current
struct inabstruct {

	double Gnab;
	double junc,sl,na;

} inab;

// Ca Background Current
struct icabstruct {

	double Gcab;
	double junc,sl,ca;

} icab;

// Background Cl current (ICaCl)
struct iclbstruct {

	double G,cftr,cl;

} iclb;

// I_K,ACh
struct ikachstruct{

	double ik;
	double fKAch_Na,RK_Ach;
	double Gmax,g;

} ikach;

// Calcium leak via SERCA pump
struct jleakstruct {

	double vleak,ca;

} jleak;

// Calcium uptake via SERCA pump
struct jupstruct {

	double vmaxup,kup,ca,rategup;

} jup;

// SR calcium release flux, via RyR (Jrel)
struct jrelstruct {
	
	double maxsr,minsr;
	double kcasr,kosrca,koca,kisrca;
	double ri,o,kim,kom,kica;
	double EC,ks,Q10SRCaP,Vmax_SRCaP,Kmf,Kmr,hillSRCaP;
	double SRCarel,SERCA,SRleak;

} jrel;

// Ixfer 
struct jxferstruct {

	double vxfer,ca;

} jxfer;

// Ca buffer
struct bufstruct {

	double kon_na,koff_na;
	double Bmax_Naj,Bmax_Nasl;
	double kon_tncl,koff_tncl;
	double Bmax_TnClow;
	double kon_tnchca,koff_tnchca;
	double Bmax_TnChigh;
	double kon_tnchmg,koff_tnchmg;
	double kon_cam,koff_cam;
	double Bmax_CaM;
	double kon_myoca,koff_myoca,kon_myomg,koff_myomg;
	double Bmax_myosin;
	double kon_sr,koff_sr;
	double Bmax_SR;
	double kon_sll,koff_sll;
	double Bmax_SLlowj,Bmax_SLlowsl;
	double kon_slh,koff_slh;
	double Bmax_SLhighj,Bmax_SLhighsl;
	double kon_csqn,koff_csqn;
	double Bmax_Csqn;

	double J_CaB_cytosol,J_CaB_junction,J_CaB_sl;
	//double cai,casr,cass;

} buf;

// Translocation of Ca Ions from NSR to JSR
struct jtrstruct {

	double tau,ca;

} jtr;

void val_consts(FILE *);
void make_ExPTable();

//void eular(int n,double h,double x[],double t);
void runge(int n,double h,double x[],double t);
void function(double x[],double f[],double t);
void input_para(FILE *);
//void initial_mem(int tMAX);
void initial_mem();
void closed_mem();

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

void data_out(FILE *, double t, double u[]);
void current(FILE *,FILE *,FILE *,FILE *,FILE *,FILE *,FILE *, double t,double x[]);

void out_ikr (FILE *, double time, double p[]);
void out_iks (FILE *, double time, double p[]);
void out_ical(FILE *, double time, double p[]);
void out_inaca (FILE *f, double time, double p[]);
void out_inak (FILE *f, double time, double p[]);
void out_cicr (FILE *f, double time, double p[]);
void out_ikach (FILE *, double time, double p[]);

void comp_reversal_potential(double x[]);
void comp_ina(double x[]);
void comp_inal(double x[]);
void comp_ito(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iks(double x[]);
void comp_ik1(double x[]);
void comp_inaca(double x[]);
void comp_inak(double x[]);
void comp_ipca(double x[]);
void comp_ikp(double x[]);
void comp_icab(double x[]);
void comp_inab(double x[]);
void comp_iclb(double x[]);
void comp_iclca(double x[]);
void comp_buffer(double x[]);
void comp_jrel(double x[]);
void comp_concentration (double x[]);
void comp_ikach(double x[]);

void current_ikr(FILE *, double t, double x[]);
void current_iks(FILE *, double t, double x[]);
void current_ical(FILE *, double t, double x[]);
void current_incx(FILE *, double t, double x[]);
void current_inak(FILE *, double t, double x[]);
void current_ik1(FILE *, double t, double x[]);
void current_ina(FILE *, double t, double x[]);
void current_ito(FILE *, double t, double x[]);
void current_it(FILE *, double t, double x[]);
void current_irel(FILE *, double t, double x[]);
void current_ikach(FILE *, double t, double x[]);

main(int argc,char **argv);
