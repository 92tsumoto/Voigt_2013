#include "syspara.h"

void comp_ikach(double x[])
{
	
	double f_K_Ach_Na,RK_Ach;
	
	f_K_Ach_Na = 1.6863/(1.0+pow(10.0/x[36],3));
	RK_Ach = 0.055 + 0.40/(1.0+exp((x[0]-var.Ek + 9.53)/17.18));

	if(var.celltype==0 && var.simtype==2){
		ikach.ik = ikach.g*(1.0+f_K_Ach_Na)*(var.CCh/(var.CCh + 0.125))*RK_Ach*(x[0]-var.Ek);
	} else if(var.celltype==1 && var.simtype==2){
		ikach.ik = ikach.g*(var.CCh/(var.CCh + 0.125))*RK_Ach*(x[0]-var.Ek);
	} else {
		ikach.ik = 0.0;
	}
	//printf("ikach=%lf\n",ikach.ik);

}

void comp_ina(double x[])
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	//printf("iV=%d,V1=%f,V2=%f,d1=%f,d2=%f\n",iV,V1,V2,d1,d2);

	ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
	ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
	ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
	ina.tauh = ina.Ttauh[iV]*d2 + ina.Ttauh[iV+1]*d1;
	ina.jss = ina.Tjss[iV]*d2 + ina.Tjss[iV+1]*d1;
	ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;

	if(ina.mss < 0.0 ) { ina.mss = 0.0; }
	if(ina.mss > 1.0 ) { ina.mss = 1.0; }

	//if(x[0] < -95.0){ x[1]=0.0;}
	//if(x[0] < -120.0){ x[2]=1.0;}
	//if(x[0] < -220.0){ x[3]=1.0;}

	ina.fast_junc = var.Fjunc*ina.Gna*(x[0]-var.Ena_junc)*x[1]*x[1]*x[1]*x[2]*x[3];
	ina.fast_sl = var.Fsl*ina.Gna*(x[0]-var.Ena_sl)*x[1]*x[1]*x[1]*x[2]*x[3];
	ina.na = ina.fast_junc + ina.fast_sl;
}

void comp_inal(double x[])
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ina.mlss = ina.Tmlss[iV]*d2 + ina.Tmlss[iV+1]*d1;
	ina.tauml = ina.Ttauml[iV]*d2 + ina.Ttauml[iV+1]*d1;
	ina.hlss = ina.Thlss[iV]*d2 + ina.Thlss[iV+1]*d1;

	if(ina.mlss < 0.0 ) { ina.mlss = 0.0; }
	if(ina.mlss > 1.0 ) { ina.mlss = 1.0; }

	if(x[0] < -95.0){ x[4]=0.0;}
	if(x[0] < -120.0){ x[5]=1.0;}

	if(var.celltype==1){
		ina.late_junc = var.Fjunc*ina.Gnal*(x[0]-var.Ena_junc)*x[4]*x[4]*x[4]*x[5];
		ina.late_sl = var.Fsl*ina.Gnal*(x[0]-var.Ena_sl)*x[4]*x[4]*x[4]*x[5];
		ina.late_na = ina.late_junc + ina.late_sl;
	} else {
		ina.late_junc = 0.0;ina.late_sl = 0.0;ina.late_na = 0.0;
	}
}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (double x[])
{
       
	ik1.ak1 = 1.0/(1.0+pow(x[36]/8.247,1.388))*1.0/(1.0+exp(0.2385*(x[0]-var.Ek-59.215)));
	ik1.bk1 = (0.49124*exp(0.08032*(x[0]-var.Ek+5.476))+exp(0.0618*(x[0]-var.Ek-594.31)))/(1.0+exp(-0.5143*(x[0]-var.Ek+4.753)));

	ik1.k1ss = ik1.ak1/(ik1.ak1 + ik1.bk1);

	ik1.ik = ik1.gk1*ik1.k1ss*(x[0]-var.Ek);

}

// Ito Transient Outward Current
void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.rss = ito.Trss[iV]*d2 + ito.Trss[iV+1]*d1;
	ito.taur = ito.Ttaur[iV]*d2 + ito.Ttaur[iV+1]*d1;

	ito.sss = ito.Tsss[iV]*d2 + ito.Tsss[iV+1]*d1;
	ito.taus = ito.Ttaus[iV]*d2 + ito.Ttaus[iV+1]*d1;

	ito.fast = ito.Gtof*x[8]*x[9]*(x[0]-var.Ek);
	ito.slow = 0.0;

	ito.ik = ito.fast+ito.slow;
}

// Ultra Rapid Activating Potassium Current 
void comp_ikur (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	// activation
	ikur.xkurss = ikur.Txkurss[iV]*d2 + ikur.Txkurss[iV+1]*d1;
	ikur.tauxkur = ikur.Ttauxkur[iV]*d2 + ikur.Ttauxkur[iV+1]*d1;
	// inactivation
	ikur.ykurss = ikur.Tykurss[iV]*d2 + ikur.Tykurss[iV+1]*d1;
	ikur.tauykur = ikur.Ttauykur[iV]*d2 + ikur.Ttauykur[iV+1]*d1;

	ikur.ik = ikur.Gkur*x[10]*x[11]*(x[0]-var.Ek);
	//printf("GKur=%lf\n",ikur.Gkur);

}

// Rapidly Activating Potassium Current 
void comp_ikr (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikr.xrss = ikr.Txrss[iV]*d2 + ikr.Txrss[iV+1]*d1;
	ikr.tauxr = ikr.Ttauxr[iV]*d2 + ikr.Ttauxr[iV+1]*d1;
	ikr.rkr = ikr.Trkr[iV]*d2 + ikr.Trkr[iV+1]*d1;

	ikr.ik = ikr.gkr*ikr.rkr*x[6]*(x[0]-var.Ek);

}

// Slowly Activating Potassium Current 
void comp_iks (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.xsss = iks.Txsss[iV]*d2 + iks.Txsss[iV+1]*d1;
	iks.tauxs = iks.Ttauxs[iV]*d2 + iks.Ttauxs[iV+1]*d1;

	iks.junc = var.Fjunc*iks.gks_junc*x[7]*x[7]*(x[0]-var.Eks);
	iks.sl = var.Fsl*iks.gks_sl*x[7]*x[7]*(x[0]-var.Eks);
	
	iks.ik = iks.junc + iks.sl;
}

// Plateu Potassium Current 
void comp_ikp (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikp.ss = ikp.Tss[iV]*d2 + ikp.Tss[iV+1]*d1;

	ikp.junc = var.Fjunc*ikp.Gkp*ikp.ss*(x[0]-var.Ek);
	ikp.sl = var.Fsl*ikp.Gkp*ikp.ss*(x[0]-var.Ek);

	ikp.ik = ikp.junc + ikp.sl;

}

// Ca-activated Cl Current 
void comp_iclca (double x[])
{
	
	iclca.junc = var.Fjunc*iclca.Gclca/(1.0+iclca.kd_cl_ca/x[38])*(x[0]-var.Ecl);
	iclca.sl = var.Fsl*iclca.Gclca/(1.0+iclca.kd_cl_ca/x[39])*(x[0]-var.Ecl);

	iclca.cl = iclca.junc + iclca.sl;

}


// L-type calcium current
void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	// VDA
	ical.dss = ical.Tdss[iV]*d2 + ical.Tdss[iV+1]*d1;
	ical.taud = ical.dss*(ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1);
	// VDI 
	ical.fss = ical.Tfss[iV]*d2 + ical.Tfss[iV+1]*d1;
	ical.tauf = ical.Ttauf[iV]*d2 + ical.Ttauf[iV+1]*d1;

	// CDI 
	//ical.fcaCaMSL = 0.1/(1.0+(0.01/x[39]));
	//ical.fcaCaj = 0.1/(1.0+(0.01/x[38]));
	ical.fcaCaMSL = 0.0;
	ical.fcaCaj = 0.0;

	// temporary val
	ical.tmp1 = ical.Ttmp1[iV]*d2 + ical.Ttmp1[iV+1]*d1;
	ical.tmp2 = ical.Ttmp2[iV]*d2 + ical.Ttmp2[iV+1]*d1;


	ical.bar_ca_j =ical.pca*4.0*(x[0]*F/var.RTonF)*(0.341*x[38]*ical.tmp1-0.341*var.cao)/(ical.tmp1-1.0);
	ical.bar_ca_sl =ical.pca*4.0*(x[0]*F/var.RTonF)*(0.341*x[39]*ical.tmp1-0.341*var.cao)/(ical.tmp1-1.0);

	ical.bar_k =ical.pk*(x[0]*F/var.RTonF)*(0.75*x[37]*ical.tmp2-0.75*var.ko)/(ical.tmp2-1.0);
	
	ical.bar_na_j =ical.pna*(x[0]*F/var.RTonF)*(0.75*x[34]*ical.tmp2-0.75*var.nao)/(ical.tmp2-1.0);
	ical.bar_na_sl =ical.pna*(x[0]*F/var.RTonF)*(0.75*x[35]*ical.tmp2-0.75*var.nao)/(ical.tmp2-1.0);

	ical.junc = (var.Fjunc_CaL*ical.bar_ca_j*x[12]*x[13]*((1.0-x[14])+ical.fcaCaj)*pow(ical.Q10CaL,ical.Qpow))*0.45;
	ical.sl = (var.Fsl_CaL*ical.bar_ca_sl*x[12]*x[13]*((1.0-x[15])+ical.fcaCaMSL)*pow(ical.Q10CaL,ical.Qpow))*0.45;
	ical.ca = ical.junc + ical.sl;

	ical.k = (ical.bar_k*x[12]*x[13]*(var.Fjunc_CaL*((1.0-x[14])+ical.fcaCaj)+var.Fsl_CaL*((1.0-x[15])+ical.fcaCaMSL))*pow(ical.Q10CaL,ical.Qpow))*0.45;

	ical.cana_j = (var.Fjunc_CaL*ical.bar_na_j*x[12]*x[13]*((1.0-x[14])+ical.fcaCaj)*pow(ical.Q10CaL,ical.Qpow))*0.45;
	ical.cana_sl = (var.Fsl_CaL*ical.bar_na_sl*x[12]*x[13]*((1.0-x[15])+ical.fcaCaMSL)*pow(ical.Q10CaL,ical.Qpow))*0.45;
	ical.cana = ical.cana_j + ical.cana_sl;

	ical.total = ical.ca + ical.k + ical.cana;

}

// Na-Ca Exchanger NCX
void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ncx.hca=ncx.Thca[iV]*d2 + ncx.Thca[iV+1]*d1;
	ncx.hna=ncx.Thna[iV]*d2 + ncx.Thna[iV+1]*d1;

	ncx.s1_junc = ncx.hca*x[34]*x[34]*x[34]*var.cao;
	ncx.s1_sl = ncx.hca*x[35]*x[35]*x[35]*var.cao;
	ncx.s2_junc = ncx.hna*var.nao*var.nao*var.nao*x[38];
	ncx.s2_sl = ncx.hna*var.nao*var.nao*var.nao*x[39];
	ncx.s3_junc = ncx.kmcai*var.nao*var.nao*var.nao*(1.0+(x[34]/ncx.kmnai)*(x[34]/ncx.kmnai)*(x[34]/ncx.kmnai))
					+ncx.kmnao*ncx.kmnao*ncx.kmnao*x[38]*(1.0+x[38]/ncx.kmcai)+ncx.kmcao*x[34]*x[34]*x[34]
					+x[34]*x[34]*x[34]*var.cao + var.nao*var.nao*var.nao*x[38];
	ncx.s3_sl = ncx.kmcai*var.nao*var.nao*var.nao*(1.0+(x[35]/ncx.kmnai)*(x[35]/ncx.kmnai)*(x[35]/ncx.kmnai))
					+ncx.kmnao*ncx.kmnao*ncx.kmnao*x[39]*(1.0+x[39]/ncx.kmcai)+ncx.kmcao*x[35]*x[35]*x[35]
					+x[35]*x[35]*x[35]*var.cao + var.nao*var.nao*var.nao*x[39];

	ncx.Ka_junc = 1.0/(1.0+(ncx.kdact/x[38])*(ncx.kdact/x[38]));
	ncx.Ka_sl = 1.0/(1.0+(ncx.kdact/x[39])*(ncx.kdact/x[39]));

	ncx.junc = var.Fjunc*ncx.Gmax*pow(ncx.Q10NCX,ical.Qpow)*ncx.Ka_junc*(ncx.s1_junc-ncx.s2_junc)/ncx.s3_junc/(1.0+ncx.ksat*ncx.hna);
	ncx.sl = var.Fsl*ncx.Gmax*pow(ncx.Q10NCX,ical.Qpow)*ncx.Ka_sl*(ncx.s1_sl-ncx.s2_sl)/ncx.s3_sl/(1.0+ncx.ksat*ncx.hna);
	ncx.j = ncx.junc + ncx.sl;

}

// Na-K Pump
void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.knai = inak.Tknai[iV]*d2 + inak.Tknai[iV+1]*d1;
	inak.knao = inak.Tknao[iV]*d2 + inak.Tknao[iV+1]*d1;

	inak.fnak = 1.0/(1.0 + inak.knai + inak.knao);

	inak.junc = var.Fjunc*inak.Gmax*inak.fnak*var.ko/(1.0+pow(inak.kmna/x[34],4.0))/(var.ko + inak.kmk);
	inak.sl = var.Fsl*inak.Gmax*inak.fnak*var.ko/(1.0+pow(inak.kmna/x[35],4.0))/(var.ko + inak.kmk);
	inak.na = inak.junc + inak.sl;

}

// Sarcolemmal Ca Pump 
void comp_ipca (double x[])
{
	ipca.junc = var.Fjunc*pow(ipca.Q10SLCaP,ical.Qpow)*ipca.Pmax*pow(x[38],1.6)/(pow(ipca.kmPca,1.6) + pow(x[38],1.6));
	ipca.sl = var.Fsl*pow(ipca.Q10SLCaP,ical.Qpow)*ipca.Pmax*pow(x[39],1.6)/(pow(ipca.kmPca,1.6) + pow(x[39],1.6));
	ipca.ca = ipca.junc + ipca.sl;
}

// Ca Background Current 
void comp_icab (double x[])
{
	icab.junc = var.Fjunc*icab.Gcab*(x[0] - var.Eca_junc);
	icab.sl = var.Fsl*icab.Gcab*(x[0] - var.Eca_sl);
	icab.ca = icab.junc + icab.sl;
	//printf("Gcab = %lf\n",icab.Gcab);
}

// Na Background Current 
void comp_inab (double x[])
{
	inab.junc = var.Fjunc*inab.Gnab*(x[0] - var.Ena_junc);
	inab.sl = var.Fsl*inab.Gnab*(x[0] - var.Ena_sl);
	inab.na = inab.junc + inab.sl;
}

// Cl Background Current 
void comp_iclb (double x[])
{
	iclb.cl = iclb.G*(x[0] - var.Ecl) + iclb.cftr*(x[0] - var.Ecl);
}

void comp_jrel (double x[])
{
	
	jrel.kcasr = jrel.maxsr - (jrel.maxsr - jrel.minsr)/(1.0+pow(jrel.EC/x[33],2.5));
	jrel.kosrca = jrel.koca/jrel.kcasr;
	jrel.kisrca = jrel.kica*jrel.kcasr;
	jrel.ri = 1.0-x[16]-x[17]-x[18];
	//printf("maxsr=%lf minsr=%lf kcasr=%lf\n",jrel.maxsr,jrel.minsr,jrel.kcasr);
	//printf("x[16]=%lf x[17]=%lf x[18]=%lf\n",x[16],x[17],x[18]);
	//printf("kosrca=%lf kisrca=%lf ri=%lf\n",jrel.kosrca,jrel.kisrca,jrel.ri);
	
	jrel.SRCarel = jrel.ks*x[17]*(x[33]-x[38]);
	jrel.SERCA = pow(jrel.Q10SRCaP,ical.Qpow)*jrel.Vmax_SRCaP*(pow(x[40]/jrel.Kmf,jrel.hillSRCaP)-pow(x[33]/jrel.Kmr,jrel.hillSRCaP))
				/(1.0+pow(x[40]/jrel.Kmf,jrel.hillSRCaP)+pow(x[33]/jrel.Kmr,jrel.hillSRCaP));
	if(var.celltype == 0){  // control (mM/ms)
		jrel.SRleak = 5.348E-6*(x[33]-x[38]);
	} else if(var.celltype == 1){	// AF case (mM/ms)
		jrel.SRleak = (1.0+0.25)*5.348E-6*(x[33]-x[38]);
	}

}

void comp_buffer (double x[])
{
	//J_CaB_cytosol = f[21]+f[22]+f[23]+f[24]+f[25]+f[26]+f[27];
	buf.J_CaB_cytosol = (buf.kon_tncl*x[40]*(buf.Bmax_TnClow-x[21])-buf.koff_tncl*x[21])
						+ (buf.kon_tnchca*x[40]*(buf.Bmax_TnChigh-x[22]-x[23])-buf.koff_tnchca*x[22])
						+ (buf.kon_tnchmg*var.mgi*(buf.Bmax_TnChigh-x[22]-x[23])-buf.koff_tnchmg*x[23])
						+ (buf.kon_cam*x[40]*(buf.Bmax_CaM-x[24])-buf.koff_cam*x[24]);
						+ (buf.kon_myoca*x[40]*(buf.Bmax_myosin-x[25]-x[26])-buf.koff_myoca*x[25])
						+ (buf.kon_myomg*var.mgi*(buf.Bmax_myosin-x[25]-x[26])-buf.koff_myomg*x[26])
						+ (buf.kon_sr*x[40]*(buf.Bmax_SR-x[27])-buf.koff_sr*x[27]);
	//J_CaB_junction = f[28]+f[30]
	buf.J_CaB_junction = (buf.kon_sll*x[38]*(buf.Bmax_SLlowj-x[28])-buf.koff_sll*x[28])
						+ (buf.kon_slh*x[38]*(buf.Bmax_SLhighj-x[30])-buf.koff_slh*x[30]);
	//J_CaB_sl = f[29]+f[31]
	buf.J_CaB_sl = (buf.kon_sll*x[39]*(buf.Bmax_SLlowsl-x[29])-buf.koff_sll*x[29]) 
					+ (buf.kon_slh*x[39]*(buf.Bmax_SLhighsl-x[31])-buf.koff_slh*x[31]);
	//printf("J_CaB_sytosol=%lf J_CaB_junction=%lf J_CaB_sl=%lf\n",buf.J_CaB_cytosol,buf.J_CaB_junction,buf.J_CaB_sl);

}

// Reversal potentials */

void comp_reversal_potential(double x[])
{
	var.Ena_junc = var.RTonF*log(var.nao/x[34]);	// [Na]_junc --> x[34]; junctional Na concentration
	var.Ena_sl = var.RTonF*log(var.nao/x[35]);		// [Na]_sl --> x[35]; salcolemmal Na concentration
	var.Ek = var.RTonF*log(var.ko/x[37]);			// [K]i --> x[37]; salcolemmal K concentration
	var.Eks = var.RTonF*log((var.ko+var.prnak*var.nao)/(x[37]+var.prnak*x[36])); // [K]i --> x[37]; [Na]i --> x[36]
	var.Eca_junc = var.RTon2F*log(var.cao/x[38]);
	var.Eca_sl = var.RTon2F*log(var.cao/x[39]);
	
	//printf("Ena=%lf, Ek=%lf, Eks=%lf\n",var.Ena,var.Ek,var.Eks);
}

