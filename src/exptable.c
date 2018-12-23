#include "syspara.h"

void make_ExpTable()
{

	int vindex,kiindex;
	double v,ki;
	double am,bm,ah,bh,ad,bd,aj,bj,aml,bml;
	double ua,ub,ia,ib;
	double axr,bxr;
	double axs,bxs;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-(double)Emax;

		// fast INa
		//if(var.model_type==0){	
		//  Grandi et al. Circ Res. 2011 model of INa
		/*
			ina.Tmss[vindex] = 1.0/((1.0+exp((-v-56.86)/9.03))*(1.0+exp((-v-56.86)/9.03)));
			ina.Ttaum[vindex] = 0.1292*exp(-((v+45.79)/15.54)*((v+45.79)/15.54)) + 0.06487*exp(-((v-4.823)/51.12)*((v-4.823)/51.12));
			
			ina.Thss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
			
			if(v<-40.0){
				ah = 0.057*exp( -(v+80.0) / 6.8 );
			} else {
				ah = 0.0;
			}
			if(v<-40.0){
				bh = 2.7*exp(0.079*v)+3.1E+5*exp(0.3485*v);
			} else {
				bh = 0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
			}
			ina.Ttauh[vindex] = 1.0/(ah+bh);
		
			ina.Tjss[vindex] = ina.Thss[vindex];	// 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
			
			if(v<-40.0){
				aj = ((-2.5428E+4*exp(0.2444*v)-6.948E-6*exp(-0.04391*v))*(v+37.78))/(1.0+exp(0.311*(v+79.23)));
			} else {
				aj = 0.0;
			}
			if(v<-40.0){
				bj = (0.02424*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
			} else {
				bj = (0.6*exp(0.057*v))/(1.0+exp(-0.1*(v+32.0)));
			}
			ina.Ttauj[vindex] = 1.0/(aj+bj);

	//	} else if(var.model_type == 1){	// Courtemanche et al. model
	*/		
			if(v==47.13){ 
				am = 3.2;
			} else {
				am = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
			}
			bm = 0.08*exp(-v/11.0);
			ina.Tmss[vindex] = am/(am+bm);
			ina.Ttaum[vindex] = 1.0/(am+bm);
			
			if(v<-40.0){
				ah = (0.135*exp(-(v+80.0)/6.8));
			} else {
				ah = 0.0;
			}
			if(v < -40){
				bh = 3.56*exp(0.079*v) + 3.1E+5*exp(0.35*v);
			} else {
				bh = (1.0/(0.13*(1.0+exp(-(v+10.66)/11.1))));
			}
			ina.Thss[vindex] = ah/(ah+bh);
			ina.Ttauh[vindex] = 1.0/(ah+bh);

			if(v<-40.0){
				aj = (((-127140*exp(0.2444*v) -3.474E-5*exp(-0.04391*v))*(v+37.78))/(1.0+exp(0.311*(v+79.23))));
			} else {
				aj = 0.0;
			}
			
			if(v<-40.0){
				bj = ((0.1212*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14))));
			} else {
				bj = ((0.3*exp(-2.5535E-7*v))/(1.0+exp(-0.1*(v+32.0))));
			}
			ina.Tjss[vindex] = aj/(aj+bj);
			ina.Ttauj[vindex] = 1.0/(aj+bj);

	//	}

		// Late INa
		if(fabs(v+47.13)<0.001){
			aml = 0.16*v+10.7408;
		} else {
			aml = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
		}
		bml = 0.08*exp(-v/11.0);
		ina.Tmlss[vindex] = aml/(aml+bml);
		ina.Ttauml[vindex] = 1.0/(aml+bml);
		ina.Thlss[vindex] = 1.0/(1.0+exp((v+91.0)/6.1));

		// for ikr 
		ikr.Txrss[vindex] = 1.0/(1.0+exp(-(v+10.0)/5.0));
		ikr.Ttauxr[vindex] = 550.0/(1.0+exp((-22.0-v)/9.0))*6.0/(1.0+exp((v+11.0)/9.0)) + 230.0/(1.0+exp((v+40.0)/20.0));
		ikr.Trkr[vindex] = 1.0/(1.0+exp((v+74.0)/24.0));

		// for iks 
		if(var.simtype==0 || var.simtype==2){	// control
			iks.Txsss[vindex] = 1.0/(1.0+exp((-3.8-v)/14.25));
			iks.Ttauxs[vindex] = 990.1/(1.0+exp(-(v+2.436)/14.12));
		} else if(var.simtype == 1){ // with ISO stimulation
			iks.Txsss[vindex] = 1.0/(1.0+exp((-40.0-3.8-v)/14.25));
			iks.Ttauxs[vindex] = 990.1/(1.0+exp(-(v+2.436+40.0)/14.12));
		}

		// for ikp 
		ikp.Tss[vindex] = 1.0/(1.0+exp((7.488-v)/5.98));

		// ito
		ito.Trss[vindex] = 1.0/(1.0+exp(-(v+1.0)/11.0));
		ito.Ttaur[vindex] = 3.5*exp(-(v/30.0)*(v/30.0))+1.5;
		
		ito.Tsss[vindex] = 1.0/(1.0+exp((v+40.5)/11.5));
		ito.Ttaus[vindex] = 25.635*exp(-((v+52.45)/15.8827)*((v+52.45)/15.8827))+24.14;

		// for ikur 
		ikur.Txkurss[vindex] = 1.0/(1.0+exp(-(v+6.0)/8.6));
		ikur.Ttauxkur[vindex] = 9.0/(1.0+exp((v+5.0)/12.0)) + 0.5;
		ikur.Tykurss[vindex] = 1.0/(1.0+exp((7.5+v)/10.0));
		ikur.Ttauykur[vindex] = 590.0/(1.0+exp((v+60.0)/10.0)) + 3050.0;

		// for ical
		if(var.simtype==0 || var.simtype==2){	// control;
			ical.Tdss[vindex] = 1.0/(1.0+exp(-(v+9.0)/6.0));
			if(fabs(v+9)<0.001){
				ical.Ttaud[vindex] = 1.19;
			} else {
				ical.Ttaud[vindex] = (1.0-exp(-(v+9.0)/6.0))/(0.035*(v+9.0));
			}
		} else if(var.simtype == 1){	// with ISO stimulation
			ical.Tdss[vindex] = 1.0/(1.0+exp(-(v+9.0+3.0)/6.0));
			if(fabs(v+9.0+3.0)<0.001){
				ical.Ttaud[vindex] = -0.397*v;	// -(25*v/63)
			} else {
				ical.Ttaud[vindex] = (1.0-exp(-(v+9.0+3.0)/6.0))/(0.035*(v+9.0+3.0));
			}
		}
		
		if(var.simtype == 0 || var.simtype==2){	// control
			ical.Tfss[vindex] = 1.0/(1.0+exp((v+30.0)/7.0))+0.2/(1.0+exp((50.0-v)/20.0));
			ical.Ttauf[vindex] = 1.0/(0.0197*exp(-(0.0337*(v+25.0))*(0.0337*(v+25.0)))+0.02);
		} else if(var.simtype == 1){
			ical.Tfss[vindex] = 1.0/(1.0+exp((v+30.0+3.0)/7.0))+0.2/(1.0+exp((50.0-v-3.0)/20.0));
			ical.Ttauf[vindex] = 1.0/(0.0197*exp(-(0.0337*(v+25.0+3.0))*(0.0337*(v+25.0+3.0)))+0.02);
		}
		ical.Ttmp1[vindex] = exp(v/var.RTon2F);
		ical.Ttmp2[vindex] = exp(v/var.RTonF);

		// inak 
		inak.Tknai[vindex] = 0.1245*exp((-0.1*v)/var.RTonF);
		inak.Tknao[vindex] = 0.0363*inak.sigma*exp((-1.0*v)/var.RTonF);

		// incx
		ncx.Thca[vindex] = exp(ncx.gamma*v/var.RTonF);
		ncx.Thna[vindex] = exp((ncx.gamma-1.0)*v/var.RTonF);
		
	} //for i loop end



}
