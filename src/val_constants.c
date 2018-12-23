
#include "syspara.h"

void val_consts(FILE *fp1)
{
	int i,w;
	double v_old,dvdt,dvdt_new;

	// Cell Geometry */
	var.length = 100;	// Length of the cell (um)
	var.a = 10.25;     	// Radius of the cell (um)
	var.vcell = 1.0E-15*M_PI*var.a*var.a*var.length; // Cell Volume (L)
	var.vmyo = var.vcell*0.65;	// Cytoplasmic volume (L)
	var.vsr  = var.vcell*0.035;	// Sarcoplasmic volume (L)
	var.vsl  = var.vcell*0.02;	// Sarcomere volume (L)
	var.vjunc  = var.vcell*0.0539*0.01;	//  junctional volume (L)
	//var.vjunc  = var.vcell*0.001;	//  junctional volume (L)

	var.Fjunc = 0.11;
	var.Fsl = 1.0-var.Fjunc;
	var.Fjunc_CaL = 0.9;
	var.Fsl_CaL = 1.0-var.Fjunc_CaL;

	var.ageo = 2.0*M_PI*var.a*var.a + 2.0*M_PI*var.a*var.length;  // eometric membrane area (um^2)
	var.acap = var.ageo*2.0;		// Capacitive membrane area --> 142.0 (pF)
	var.Cmem = 1.1E-10;			// (F) membrane capacitance 1.3810E-10 --> 138.1 (pF)

	var.J_na_juncsl = 1.0/(1.6382E+12/3.0*100.0);	// (L/ms)=6.1043E-13
	var.J_na_slmyo = 1.0/(1.8308E+10/3.0*100.0);	// (L/ms)=5.4621E-11
	var.J_ca_juncsl = 1.0/1.2134E+12;				// (L/ms)=8.2413e-13
	var.J_ca_slmyo = 1.0/2.68510E+11;				// (L/ms)=3.2743E-12 (miss?) <-- 3.7243E-12

	// Q10
	//var.K_Q10 = 3.0;

	// Ion Valences
	var.zna = 1.0;  // Na valence
	var.zk = 1.0;   // K valence
	var.zca = 2.0;  // Ca valence

	// invariant constant
	var.RTonF = R*T/F;
	var.RTon2F = R*T/(var.zca*F);
	var.RTon4F = R*T/(var.zca*var.zca*F*F);

	// Extracellular Concentrations
	var.nao = 140.0;     // Initial Bulk Medium Na (mM)
	var.ko = 5.4;      // Initial Bulk Medium K (mM)
	var.cao = 1.8;     // Initial Bulk Medium Ca (mM)
	var.clo = 150.0;     // Initial Bulk Medium Cl (mM)
	var.mgi = 1.0;     // Initial intracellular Mg (mM)
	var.cli = 15.0;     // Initial intracellular Cl (mM)
	var.Ecl = var.RTonF*log(var.cli/var.clo);

	// Fast sodium current
		if(var.celltype == 0){	// control
			//ina.Gna = 23.0;	// (nS/pF).
			ina.Gna = 7.8;	// (nS/pF).
		} else if(var.celltype == 1){	// AF case
			//ina.Gna = 0.9*23.0;	// (nS/pF).
			ina.Gna = 0.9*7.8;	// (nS/pF).
		}
		
	// Late sodium current
		if(var.celltype == 0){ 	// control
			ina.Gnal = 0.0025*0.0;	// (nS/pF).
		} else if(var.celltype == 1){	// AF case
			ina.Gnal = 0.0025*1.0;	// (nS/pF).
		}
		ina.tauhl = 600.0;
		
	// Inward rectifier K current: Ik1
		if(var.celltype==0){
			ik1.Gk1max = 0.0525; // (nS/pF)
		} else if(var.celltype==1){ 	// AF
			ik1.Gk1max = 0.0525*1.8; // (nS/pF)
		}
		ik1.gk1 = ik1.Gk1max*sqrt(var.ko/5.4);

	// Transient outward current
		if(var.celltype==0){	//control
			ito.Gtof = 0.165;	// (nS/pF).
		} else if(var.celltype==1){	// AF case 
			ito.Gtof = 0.3*0.165;	// (nS/pF).
		}

	// Ultra Rapid delayed rectifier potassium current (Ikur)
		if(var.celltype==0 && (var.simtype==0 || var.simtype==2)){	// control
			ikur.Gkur = 0.045;  //(nS/pF)
		} else if(var.celltype==1 && (var.simtype==0 || var.simtype==2)){	// AF 
			ikur.Gkur = 0.045*0.5;  //(nS/pF)
		} else if(var.celltype==2 && (var.simtype==0 || var.simtype==2)){	// Right Atrirum 
			ikur.Gkur = 0.045*1.2;  //(nS/pF)
		} else if(var.celltype==0 && var.simtype==1){	// with ISO stimulation 
			ikur.Gkur = 0.045*3.0;  //(nS/pF)
		} else if(var.celltype==1 && var.simtype==1){	// AF with ISO stimulation 
			ikur.Gkur = 0.045*0.5*3.0;  //(nS/pF)
		} else if(var.celltype==2 && var.simtype==1){	// AF with ISO stimulation 
			ikur.Gkur = 0.045*0.5*3.0*1.2;  //(nS/pF)
		}

	// Rapid delayed rectifier potassium current (Ikr)
		ikr.Gkr = 0.035;  //(nS/pF)
		ikr.gkr = ikr.Gkr*sqrt(var.ko/5.4);

	// Slow delayed rectifier potassium current (Iks)
		if(var.celltype==0 && (var.simtype==0 || var.simtype==2)){	//control
			iks.Gks = 0.0035;
		} else if(var.celltype==0 && var.simtype==1){	// with ISO stimulation 
			iks.Gks = 3.0*0.0035; 
		} else if(var.celltype==1 && (var.simtype==0 || var.simtype==2)){	// AF case
			iks.Gks = 2.0*0.0035; 
		} else if(var.celltype==1 && var.simtype==1){	// AF case with ISO stimulation
			iks.Gks = 4.0*0.0035; 
		}
		iks.gks_junc = iks.Gks; 
		iks.gks_sl = iks.Gks; 
		var.prnak = 0.01833;
	
	// ACh activated potassium current (Ik,ACh)
		if(var.celltype==0 && (var.simtype==0 || var.simtype==2)){	//control
			ikach.Gmax = 0.1275;
		} else if(var.celltype==1 && (var.simtype==0 || var.simtype==2)){	// AF case
			ikach.Gmax = 0.0691; 
		}
		ikach.g = ikach.Gmax*sqrt(var.ko/5.4);

	// plateau K Current 
		ikp.Gkp = 0.002;		// Max. conductance of plateau K current(nS/pF)

	// Ca-activated Cl current
		iclca.Gclca = 0.0548;	// (mS/uF)
		iclca.kd_cl_ca = 100E-3;	// (mM)

	// L-type calcium current
		if(var.celltype == 0 && (var.simtype==0 || var.simtype==2)){ // control
			ical.pca = 2.7E-4;	// (cm/s)
			ical.pk = 1.35E-7;	// (cm/s)
			ical.pna = 0.75E-8;	// (cm/s)
		} else if(var.celltype == 0 && var.simtype == 1){ // with ISO stimuli
			ical.pca = 2.7E-4*1.5;	// (cm/s)
			ical.pk = 1.35E-7*1.5;	// (cm/s)
			ical.pna = 0.75E-8*1.5;	// (cm/s)
		} else if(var.celltype == 1 && (var.simtype == 0 || var.simtype==2)){ // AF case
			ical.pca = 2.7E-4*0.5;	// (cm/s)
			ical.pk = 1.35E-7*0.5;	// (cm/s)
			ical.pna = 0.75E-8*0.5;	// (cm/s)
		} else if(var.celltype == 1 && var.simtype == 1){ // AF case with ISO stimuli
			ical.pca = 2.7E-4*1.5*0.5;	// (cm/s)
			ical.pk = 1.35E-7*1.5*0.5;	// (cm/s)
			ical.pna = 0.75E-8*1.5*0.5;	// (cm/s)
		}
		ical.Q10CaL = 1.8;
		ical.Qpow = (T-310.0)/10.0;

	// Na/K Pump (NaK)
	if(var.simtype==0 || var.simtype==2){	// control or Ach stimulation
		inak.kmna = 11.0;	// (mM)
	} else if(var.simtype==1){	// with ISP stimulation
		inak.kmna = 11.0*(1.0-0.25);	// (mM)
	}
		inak.kmk = 1.5;		// (mM)
		inak.Gmax = 1.26;	// (uA/uF)
		inak.Q10nak = 1.63;
		inak.Q10kmnai = 1.39;
		inak.sigma = (exp(var.nao/67.3)-1.0)/7.0;

	// Sodium-Calcium Exchanger (NCX) 
		if(var.celltype == 0){	// control
			ncx.Gmax = 3.15;	// (uA/uF)
		} else if(var.celltype == 1){	// AF case
			ncx.Gmax = 3.15*1.4;	
		}
		ncx.kdact = 0.384E-3;
		ncx.kmcai = 3.59E-3;	// (mM)
		ncx.kmcao = 1.3;	// (mM)
		ncx.kmnai = 12.29;	// (mM)
		ncx.kmcao = 87.5;	// (mM)
		ncx.ksat = 0.27;
		ncx.gamma = 0.35;
		ncx.Q10NCX = 1.57;
		//ncx.a = 2.5;

	// Sarcolemmal Ca Pump
		ipca.Pmax = 0.0471;		// (uA/uF)
		ipca.kmPca = 0.5E-3;	// Half-saturation concentration of sarcolemmal Ca pump (mM)
		ipca.Q10SLCaP = 2.35;

	// Ca Background Current 
		icab.Gcab = 6.0643E-4; // (uA/uF)

	// Na Background Current 
		inab.Gnab = 0.597E-3; //(mS/uF)

	// Cl Background Current 
		iclb.G = 9.0E-3; //(mS/uF)
		if(var.simtype == 0 || var.simtype==2){
			iclb.cftr = 0.0;
		} else if(var.simtype == 1){
			iclb.cftr = 4.9E-3;
		}

	// calcium uptake via SERCA pump (Jup)
		jup.vmaxup = 0.6*0.006375;	// Maximal Iup conductance (mM)
		jup.kup =    0.00025;	// (mM)

	// SR calcium release flux, via RyR (Jrel)
		jrel.hillSRCaP = 1.787;	// (mM) 
		if(var.simtype==0 || var.simtype==2){	// control
			jrel.Kmf = 2.5*0.246E-3;	// (mM) 
		} else if(var.simtype == 1){	// with ISO stimulation
			jrel.Kmf = (2.5-1.25)*0.246E-3;	// (mM)
		}
		jrel.Kmr = 1.7;		// (mM)
		jrel.ks = 25.0;		// (1/ms)
		jrel.EC = 0.45;		// Ca_sr half-saturation constant of kcasr (mM)
		jrel.maxsr = 15.0;
		jrel.minsr = 1.0;
		if(var.celltype == 0 && (var.simtype == 0 || var.simtype==2)){	// control
			jrel.koca = 10.0;	// (1/mM/mM/ms)
		} else if(var.celltype == 0 && var.simtype == 1){ // with ISO stimulation
			jrel.koca = 10.0+10.0;
		} else if(var.celltype == 1 && (var.simtype == 0 || var.simtype==2)){ // AF case 
			jrel.koca = 10.0+20.0;
		} else if(var.celltype == 1 && var.simtype == 1){ // AF case with ISO stimulation
			jrel.koca = 10.0+20.0;
		}
		jrel.kica = 0.5;	// (1/mM/ms)
		jrel.kim = 0.005;	// (1/ms)
		jrel.kom = 0.06;	// (1/ms)
		jrel.Q10SRCaP = 2.6;
		jrel.Vmax_SRCaP = 5.3114E-3;	// (mM/msec)-->(286 uM/L cytosol/sec)

	// Buffering parameters
		buf.kon_na = 0.1E-3;	// (1/mM/ms)
		buf.koff_na = 1.0E-3;	// (1/ms)
		buf.Bmax_Naj = 7.561;	// (mM)
		buf.Bmax_Nasl = 1.65;	// (mM)
		// TnC low affinity
		buf.kon_tncl = 32.7;	// (1/mM/ms)
		if(var.simtype==0 || var.simtype==2){
			buf.koff_tncl = 19.6E-3;	// (1/ms)
		} else if(var.simtype==1){
			buf.koff_tncl = 1.5*19.6E-3;	// (1/ms)
		}
		buf.Bmax_TnClow = 70.0E-3;	// (mM)
		// Tnc high affinity
		buf.kon_tnchca = 2.37;	// (1/mM/ms)
		buf.koff_tnchca = 0.032E-3;	// (1/ms)
		buf.Bmax_TnChigh = 140.0E-3;	// (mM)
		buf.kon_tnchmg = 3.0E-3;	// (1/mM/ms)
		buf.koff_tnchmg = 3.33E-3;	// (1/ms)
		// CaM buffering
		buf.kon_cam	= 34.0;			// (1/mM/ms)
		buf.koff_cam = 238.0E-3;	// (1/ms)
		buf.Bmax_CaM = 24.0E-3;		// (mM)
		// Myosin buffering
		buf.kon_myoca = 13.8;
		buf.koff_myoca = 0.46E-3;
		buf.Bmax_myosin = 140.0E-3; // (mM)
		buf.kon_myomg = 0.0157;
		buf.koff_myomg = 0.057E-3;
		buf.kon_sr = 100.0;
		buf.koff_sr = 60.0E-3;
		buf.Bmax_SR = 19.0*0.9E-3;
		// Junctional and SL Ca Buffers
		buf.kon_sll = 100.0;
		buf.koff_sll = 1300.0E-3;
		buf.Bmax_SLlowj = 4.6E-3*var.vmyo/var.vjunc*0.1;
		buf.Bmax_SLlowsl = 37.4E-3*var.vmyo/var.vsl;
		buf.kon_slh = 100.0;
		buf.koff_slh = 30.0E-3;
		buf.Bmax_SLhighj = 1.65E-3*var.vmyo/var.vjunc*0.1;
		buf.Bmax_SLhighsl = 13.4E-3*var.vmyo/var.vsl;
		// Csqn buffering
		buf.kon_csqn = 100.0;
		buf.koff_csqn = 65.0;
		buf.Bmax_Csqn = 140.0E-3*var.vmyo/var.vsr;

}

