#include "syspara.h"

void function(double x[],double f[],double t)
{

	int i;
	comp_reversal_potential(x);
	comp_ina(x);
	comp_inal(x);
	comp_inab(x);
	comp_ikr(x);
	comp_iks(x);
	comp_ikp(x);
	comp_ik1(x);
	comp_ito(x);
	comp_ikur(x);
	if(var.simtype==2) comp_ikach(x);
	comp_inak(x);
	comp_iclca(x);
	comp_ical(x);
	comp_icab(x);
	comp_inaca(x);
	comp_ipca(x);
	comp_iclb(x);
	comp_jrel(x);
	comp_buffer(x);

		if(var.deb==1){ 
            printf("time=%lf Istim=%lf ",t,var.Istim);
			for(i=0;i<NN;i++){printf("%e ", x[i]);}
            printf("ENa_junc=%lf, ENa_sl=%lf, EK=%lf ECa_junc=%lf ECa_sl=%lf Ecl=%lf\n",
				var.Ena_junc,var.Ena_sl,var.Ek,var.Eca_junc,var.Eca_sl,var.Ecl);
        }

	ina.tot_junc = ina.fast_junc + ina.late_junc + inab.junc + 3.0*ncx.junc + 3.0*inak.junc + ical.cana_j;
	ina.tot_sl = ina.fast_sl + ina.late_sl + inab.sl + 3.0*ncx.sl + 3.0*inak.sl + ical.cana_sl;
	//ina.tot_junc2 = 3.0*ncx.junc + 3.0*inak.junc + ical.cana_j;
	//ina.tot_sl2 = 3.0*ncx.sl + 3.0*inak.sl + ical.cana_sl;

	if(var.simtype==2){
		ik.tot = ito.ik + ikr.ik + iks.ik + ik1.ik -2.0*inak.na + ical.k + ikp.ik + ikur.ik + ikach.ik;
	} else {
		ik.tot = ito.ik + ikr.ik + iks.ik + ik1.ik -2.0*inak.na + ical.k + ikp.ik + ikur.ik;
	}
	ica.tot_junc = ical.junc + icab.junc + ipca.junc - 2.0*ncx.junc;
	ica.tot_sl = ical.sl + icab.sl + ipca.sl - 2.0*ncx.sl;
	
	icl.tot = iclca.cl + iclb.cl;

	var.Itotal = ina.tot_junc + ina.tot_sl + ik.tot + ica.tot_junc + ica.tot_sl + icl.tot;

	f[0] = -(var.Itotal+var.Istim);
	//Fast sodium current
	f[1] = (ina.mss - x[1])/ina.taum;
	f[2] = (ina.hss - x[2])/ina.tauh;
	f[3] = (ina.jss - x[3])/ina.tauj;
	//Late sodium current
	f[4] = (ina.mlss - x[4])/ina.tauml;
	f[5] = (ina.hlss - x[5])/ina.tauhl;
	// Ikr
	f[6] = (ikr.xrss - x[6])/ikr.tauxr;
	// Iks
	f[7] = (iks.xsss - x[7])/iks.tauxs;
	// Ito,f
	f[8] = (ito.rss - x[8])/ito.taur;
	f[9] = (ito.sss - x[9])/ito.taus;
	// IKur
	f[10] = (ikur.xkurss - x[10])/ikur.tauxkur;
	f[11] = (ikur.ykurss - x[11])/ikur.tauykur;
	// LTCC
	f[12] = (ical.dss - x[12])/ical.taud;
	f[13] = (ical.fss - x[13])/ical.tauf;
	// fCa_junc
	f[14] = 1.7*x[38]*(1.0 - x[14])-11.9E-3*x[14];
	// fCa_sl
	f[15] = 1.7*x[39]*(1.0 - x[15])-11.9E-3*x[15];
	// Jrel
	f[16] = (jrel.kim*jrel.ri-jrel.kisrca*x[38]*x[16])-(jrel.kosrca*x[38]*x[38]*x[16] -jrel.kom*x[17]);	// R
	f[17] = (jrel.kosrca*x[38]*x[38]*x[16]-jrel.kom*x[17])-(jrel.kisrca*x[38]*x[17]-jrel.kim*x[18]); 	// O
	f[18] = (jrel.kisrca*x[38]*x[17] - jrel.kim*x[18])-(jrel.kom*x[18]-jrel.kosrca*x[38]*x[38]*jrel.ri);	// I
	// Na and Ca buffering
	f[19] = buf.kon_na*x[34]*(buf.Bmax_Naj-x[19])-buf.koff_na*x[19];	// NaBj	 (mM/ms)
	f[20] = buf.kon_na*x[35]*(buf.Bmax_Nasl-x[20])-buf.koff_na*x[20];	// NaBsl (mM/ms)
	// Cytosolic Ca buffers
	f[21] = buf.kon_tncl*x[40]*(buf.Bmax_TnClow-x[21])-buf.koff_tncl*x[21];	// TnCL (mM/ms)
	f[22] = buf.kon_tnchca*x[40]*(buf.Bmax_TnChigh-x[22]-x[23])-buf.koff_tnchca*x[22];	// TnCHc (mM/ms)
	f[23] = buf.kon_tnchmg*var.mgi*(buf.Bmax_TnChigh-x[22]-x[23])-buf.koff_tnchmg*x[23];	// TnCHm (mM/ms)
	f[24] = buf.kon_cam*x[40]*(buf.Bmax_CaM-x[24])-buf.koff_cam*x[24];	// CaM (mM/ms)
	f[25] = buf.kon_myoca*x[40]*(buf.Bmax_myosin-x[25]-x[26])-buf.koff_myoca*x[25];	// Myosin_Ca (mM/ms)
	f[26] = buf.kon_myomg*var.mgi*(buf.Bmax_myosin-x[25]-x[26])-buf.koff_myomg*x[26];	// Myosin_Mg (mM/ms)
	f[27] = buf.kon_sr*x[40]*(buf.Bmax_SR-x[27])-buf.koff_sr*x[27];	// SRB (mM/ms)
	// Junctional and SL Ca Buffers
	f[28] = buf.kon_sll*x[38]*(buf.Bmax_SLlowj-x[28])-buf.koff_sll*x[28];	//SLLj (mM/ms)
	f[29] = buf.kon_sll*x[39]*(buf.Bmax_SLlowsl-x[29])-buf.koff_sll*x[29];	// SLLsl (mM/ms)
	f[30] = buf.kon_slh*x[38]*(buf.Bmax_SLhighj-x[30])-buf.koff_slh*x[30];	//SLHj (mM/ms)
	f[31] = buf.kon_slh*x[39]*(buf.Bmax_SLhighsl-x[31])-buf.koff_slh*x[31];	// SLHsl (mM/ms)
	// Csqn buffering
	f[32] = buf.kon_csqn*x[33]*(buf.Bmax_Csqn-x[32])-buf.koff_csqn*x[32];	// Csqn (mM/ms)
	// [Ca]sr
	f[33] = jrel.SERCA-(jrel.SRleak*var.vmyo/var.vsr + jrel.SRCarel)-(buf.kon_csqn*x[33]*(buf.Bmax_Csqn-x[32])-buf.koff_csqn*x[32]);
	// [Na]j
	f[34] = -ina.tot_junc*var.Cmem/(var.vjunc*F) + var.J_na_juncsl/var.vjunc*(x[35]-x[34]) 
			- (buf.kon_na*x[34]*(buf.Bmax_Naj-x[19])-buf.koff_na*x[19]);	// [Na]j
	// [Na]sl
	f[35] = -ina.tot_sl*var.Cmem/(var.vsl*F) + var.J_na_juncsl/var.vsl*(x[34]-x[35]) 
			+ var.J_na_slmyo/var.vsl*(x[36]-x[35]) - (buf.kon_na*x[35]*(buf.Bmax_Nasl-x[20])-buf.koff_na*x[20]);	// [Na]sl
	// [Na]i
	f[36] = var.J_na_slmyo/var.vmyo*(x[35]-x[36]);
	// [K]i
	//f[37] = -ik.tot*var.Cmem/(var.vmyo*F);
	f[37] = 0.0;
	// [Ca]junc
	f[38] = -ica.tot_junc*var.Cmem/(var.vjunc*2.0*F) + var.J_ca_juncsl/var.vjunc*(x[39]-x[38]) - buf.J_CaB_junction 
			+ jrel.SRCarel*var.vsr/var.vjunc + jrel.SRleak*var.vmyo/var.vjunc;
	// [Ca]sl
	f[39] = -ica.tot_sl*var.Cmem/(var.vsl*2.0*F) + var.J_ca_juncsl/var.vsl*(x[38]-x[39]) + var.J_ca_slmyo/var.vsl*(x[40]-x[39]) - buf.J_CaB_sl;
	// [Ca]i
	f[40] = -jrel.SERCA*var.vsr/var.vmyo - buf.J_CaB_cytosol + var.J_ca_slmyo/var.vmyo*(x[39]-x[40]);

	//for(i=0;i<NN;i++){
	//	printf("f[%d]=%e\n",i,f[i]);
	//}
}
