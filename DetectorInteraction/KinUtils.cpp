#include <iostream>
#include <cmath>
#include "KinUtils.h"

using namespace std;

KinUtils::KinUtils(const double &m_Ebeam, const double &m_Mchi, const double &m_Maprime, const double &m_Msplit, const double &m_Epsilon, const double &m_alphaD, const double &m_Ethr, const double &m_Pthr, const double &m_Pbinding, const int &m_Seed) :
		Ebeam(m_Ebeam), Mchi(m_Mchi), Maprime(m_Maprime), Msplit(m_Msplit), Seed(m_Seed), Ethr(m_Ethr), Pthr(m_Pthr), Pbinding(m_Pbinding), Epsilon(m_Epsilon), AlphaDark(m_alphaD), Alpha(1. / 137.) {

	Rand.SetSeed(m_Seed);
}

void KinUtils::Setup() {
	f_chipXsection = new TF1*[nFunctionsElastic];
	f_chieXsection = new TF1*[nFunctionsElastic];
	f_chinuclXsection = new TF1*[nFunctionsElastic];

	/*This part is really critical. The cross-section is a function of the final state recoil energy, and as "parameter" needs the incoming
	 * chi energy. However, if event-by-event we set this parameter, and then call a "GetRandom", we trigger the integration computation for ALL events,
	 * and this is very time-consuming. Instead, if the integration is already performed, a second integration (even in a different range) is not time-consuming.
	 * Instead, I "bin" wrt the incoming chi kinetic energy (nFunctionsElastic bins), and then extract the random number from the function defining the energy in that bin!
	 *
	 * MINIMUM value of the recoil total energy: its mass + threshold + binding energy
	 * MAXIMUM value: the actual form of the maximum value is complicated (see below), but a very upper limit is: primary e- beam energy + its mass - chi mass.
	 * A.C. 13/3/2019: shouldn't it be: chi energy in this bin + chi mass - e- mass??? Not critical, since the actual integral is computed below in the xsection functions..
	 *
	 * It should be fine to use this MAX value to cache integrals, the actual max value is seen below.
	 */

	cout << "KinUtils::KinUtils setting the cross-section functions" << endl;
	fflush(stdout);
	double Trmax_protonElastic, Trmax_electronElastic,thisEbeam;
	for (int ii = 0; ii < nFunctionsElastic; ii++) {

		f_chipXsection[ii] = new TF1(Form("f_chipXsection_%i", ii), this, &KinUtils::Er_chipXsection, Pthr + Pbinding + Mn, (ii + 1) * Ebeam / nFunctionsElastic + Mn - Mchi, 1);
		f_chipXsection[ii]->SetNpx(1000);
		f_chipXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / nFunctionsElastic);

		f_chieXsection[ii] = new TF1(Form("f_chieXsection_%i", ii), this, &KinUtils::Er_chieXsection, Ethr + Me, (ii + 1) * Ebeam / nFunctionsElastic + Me - Mchi, 1);
		f_chieXsection[ii]->SetNpx(1000);
		f_chieXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / nFunctionsElastic);

		//Here the VARIABLE is the KINETIC energy of the scattered nucleus. We integrate this from threshold to Ebeam-Mchi
		f_chinuclXsection[ii] = new TF1(Form("f_chinuclXsection_%i", ii), this, &KinUtils::Er_chinuclXsection, Nuclthr, (ii + 1) * Ebeam / nFunctionsElastic - Mchi, 1);
		f_chinuclXsection[ii]->SetNpx(1000);
		f_chinuclXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / nFunctionsElastic);

		/*Cache the integrals*/
		thisEbeam= (ii + 1) * Ebeam / nFunctionsElastic;

		if ((Pbinding + Pthr) < (thisEbeam - Mchi)) f_chipXsection[ii]->Integral(Pthr + Pbinding + Mn, thisEbeam + Mn - Mchi);
		if (Ethr < (thisEbeam - Mchi)) f_chieXsection[ii]->Integral(Ethr + Me, thisEbeam + Me - Mchi);
		if (Nuclthr < (thisEbeam - Mchi)) f_chinuclXsection[ii]->Integral(Nuclthr, thisEbeam - Mchi);


	}

	cout << "KinUtils::KinUtils created" << endl;
	fflush(stdout);
}

void KinUtils::PrintParameters() {
	cout << endl;
	cout << "KinUtils Parameters:" << endl;
	cout << "e- beam (GeV): \t\t" << Ebeam << endl;
	cout << "Maprime (GeV): \t \t" << Maprime << endl;
	cout << "Mchi (GeV): \t \t" << Mchi << endl;
	cout << "Split (GeV): \t  \t" << Msplit << endl;
	cout << "Epsilon: \t \t" << Epsilon << endl;
	cout << "AlphaDark: \t \t" << AlphaDark << endl;
	cout << "e- threshold: \t \t" << Ethr << endl;
	cout << "p threshold: \t \t" << Pthr << endl;

	cout << "Nuclear mass: \t \t" << Mnucl << endl;
	cout << "Nuclear charge: \t\t" << Znucl << endl;
	cout << "Nuclear threshold: \t \t" << Nuclthr << endl;
}

/*double Tr_chipXsection
 Returns the chi p -> chi p DIFFERENTIAL cross section, dsigma/dEr, where Er is the recoil proton TOTAL energy
 Parameters: x[0]:   Recoil proton total energy
 par[0]: E0, total energy of incoming chi
 */
double KinUtils::Er_chipXsection(double *x, double *par) {
	double Er = x[0];
	double E0 = par[0];

	//The reaction is: chi(p1)+N(p2)->chi(k1)+N(k2).
	//Perform the cross-section in the lab frame, where p1=(0,0,P,E), p2=(0,0,0,Mn), k2=(xx,yy,zz,Er)
	//1: Compute Lorentz-invariant dot products

	double p1p1 = Mchi * Mchi;
	double k1k1 = Mchi * Mchi;
	double p1p2 = E0 * Mn;
	double p1k1 = -0.5 * (2 * Mn * Mn - 2 * Mchi * Mchi - 2 * Er * Mn); //A.C. Ask Eder to check
	double p1k2 = Mn * (E0 + Mn - Er);
	double p2k1 = Mn * (E0 + Mn - Er);
	double p2k2 = Er * Mn;
	double k1k2 = Mn * E0;

	double k1p2 = p2k1;
	double k2p2 = p2k2;
	double k1p1 = p1k1;

	double t = 2 * Mn * Mn - 2 * Mn * Er;
	double s = Mn * Mn + Mchi * Mchi + 2 * Mn * E0;

	//2: the amplitude squared
	double ampsq = 1.0 * (k1k2 * p1p2 + p1k2 * k1p2 - Mchi * Mchi * k2p2 - Mn * Mn * k1p1 + 2 * Mchi * Mchi * Mn * Mn) / (pow((t - Maprime * Maprime), 2));
	if (ampsq < 0.0) {
		//  cout<<"Er_chipXsection: negative amplitude!!!"<<endl;
		return 0;
	}
	//3: the 2 momenta in the CM frame (before/after). Since this is elastic scattering, they're the same!!!

	double p = ((pow(s - Mchi * Mchi - Mn * Mn, 2) - 4 * Mchi * Mchi * Mn * Mn) / (4 * s));
	p = sqrt(p);
	double k = p;

	//4: the Jacobian for the transformation cos(theta_CM) --> Recoil total energy in LAB frame
	double dcostheta = Mn / (p * k); //AC: checked this formula.

	//5: the phase-space
	double S = (k / (s * p)) * dcostheta;

	//7: the Form-Factor, using a dipole approximation.
	double FF = 1 / pow(1 + (Er * Er - Mn * Mn) / Mn * Mn, 2);

	//8: put everything together

	double dsigma = pow(FF, 2) * ampsq * S; //in "GeV^-2 units" (and no coupling yet)
	dsigma = dsigma * Epsilon * Epsilon * 4 * PI * Alpha * AlphaDark;
	dsigma *= GeVm2cm2;

	return dsigma;

}

/*double Er_chieXsection
 Returns the chi e -> chi e DIFFERENTIAL cross section, dsigma/dEr, where Er is the recoil electron total energy
 Parameters: x[0]: recoil electron total energy Er
 par[0]: E0, kinetic energy of incoming chi
 */
double KinUtils::Er_chieXsection(double *x, double *par) {
	double Er = x[0];
	double E0 = par[0];

	//The reaction is: chi(p1)+e(p2)->chi(k1)+e(k2).
	//Perform the cross-section in the lab frame, where p1=(0,0,P,E), p2=(0,0,0,Me), k2=(xx,yy,zz,Er)
	//1: Compute Lorentz-invariant dot products

	double p1p1 = Mchi * Mchi;
	double k1k1 = Mchi * Mchi;
	double p1p2 = E0 * Me;
	double p1k1 = -0.5 * (2 * Me * Me - 2 * Mchi * Mchi - 2 * Er * Me);
	double p1k2 = Me * (E0 + Me - Er);
	double p2k1 = Me * (E0 + Me - Er);
	double p2k2 = Er * Me;
	double k1k2 = Me * E0;

	double k1p2 = p2k1;
	double k2p2 = p2k2;
	double k1p1 = p1k1;

	double t = 2 * Me * Me - 2 * Me * Er;
	double s = Me * Me + Mchi * Mchi + 2 * Me * E0;

	//2: the amplitude squared
	double ampsq = 1.0 * (k1k2 * p1p2 + p1k2 * k1p2 - Mchi * Mchi * k2p2 - Me * Me * k1p1 + 2 * Mchi * Mchi * Me * Me) / (pow((t - Maprime * Maprime), 2));
	if (ampsq < 0.0) {
		//  cout<<"Er_chieXsection: negative amplitude!!! "<<E0<<" "<<Er<<endl;
		return 0;
	}
	//3: the 2 momenta in the CM frame (before/after). Since this is elastic scattering, they're the same!!!
	double p = ((pow(s - Mchi * Mchi - Me * Me, 2) - 4 * Mchi * Mchi * Me * Me) / (4 * s));
	p = sqrt(p);
	double k = p;
	//4: the Jacobian for the transformation cos(theta_CM) --> Recoil total energy in LAB frame
	double dcostheta = Me / (p * k); //AC: checked this formula.

	//5: the phase-space
	double S = (k / (s * p)) * dcostheta;

	//8: put everything together
	double dsigma = ampsq * S; //in "GeV^-2 units" (and no coupling yet)
	dsigma = dsigma * Epsilon * Epsilon * 4 * PI * Alpha * AlphaDark;
	dsigma *= GeVm2cm2;
	return dsigma;
}

/*double Er_chinuclXsection
 Returns the chi Nucl -> chi nucl DIFFERENTIAL cross section, dsigma/dTr, where Tr is the recoil nucleus KINETIC energy (note the difference wrt previous cases)
 Parameters: x[0]: recoil nucleus kinetic energy Tr
 par[0]: E0, kinetic energy of incoming chi
 */
double KinUtils::Er_chinuclXsection(double *x, double *par) {
	double Tr = x[0];
	double E0 = par[0];

	//1: Use eq. 38 of 1607.01390
	double dsigma;
	dsigma = 8 * PI * Alpha * AlphaDark * Epsilon * Epsilon * Mnucl * Znucl * Znucl;
	dsigma = dsigma / (pow(Maprime * Maprime + 2 * Mnucl * Tr, 2)); // this is in GeV^-3

	//2: cconvert in cm2 / GeV
	dsigma *= GeVm2cm2;
	return dsigma;
}

/*This is the routine that, given the cross-section for the ELASTIC interaction chi-e --> chi-e OR chi-p --> chi-p
 * estracts the recoil kinematics.
 * procID is an integer, for the process under study:
 *
 * 22/04/2016: this function now correct for the proton binding energy!
 * Returns: the total cross-section, in cm2, for the interaction process, integrated over final state (with the kinetic energy cut)
 */

double KinUtils::doElasticRecoil(const TLorentzVector &chi, TLorentzVector &recoil, TLorentzVector &recoil_chi, const int &procID) {

	TVector3 v0, v1, v2;
	TVector3 p0, pr, pchi;
	double E0, Er, Echi;
	double Tr_max; //this is the maximum recoil KINETIC energy for this incoming chi
	double P0, Pr, Pchi;
	double ctheta_r, stheta_r, phi_r, sigma;
	int ii;
	E0 = chi.E();
	P0 = chi.P();
	p0 = chi.Vect();

	/*1: extract the recoil total energy from the cross-section*/
	ii = 0;
	ii = (int) (E0 / (Ebeam / nFunctionsElastic));
	if (ii >= nFunctionsElastic) ii = nFunctionsElastic - 1; //should not happen!!!

	if (procID == Proc_Pelastic) {
		Tr_max = (2 * Mn * (E0 * E0 - Mchi * Mchi)) / (2 * E0 * Mn + Mn * Mn + Mchi * Mchi); //maximum energy transfer
		if (Tr_max < (Pthr + Pbinding)) return 0; //this event is not compatible with the threshold, it is useless to proceed further
	} else if (procID == Proc_Eelastic) {
		Tr_max = (2 * Me * (E0 * E0 - Mchi * Mchi)) / (2 * E0 * Me + Me * Me + Mchi * Mchi); //maximum energy transfer
		if (Tr_max < Ethr) {
			//cout<<"THR IS: "<<Ethr<<" "<<Tr_max<<" "<<E0<<endl;
			return 0; //this event is not compatible with the threshold, it is useless to proceed further
		}
	} else if (procID == Proc_Nuclelastic) {
		Tr_max = (2 * Mnucl * (E0 * E0 - Mchi * Mchi)) / (2 * E0 * Mnucl + Mnucl * Mnucl + Mchi * Mchi); //maximum energy transfer
		if (Tr_max < Nuclthr) {
			//cout<<"THR IS: "<<Ethr<<" "<<Tr_max<<" "<<E0<<endl;
			return 0; //this event is not compatible with the threshold, it is useless to proceed further
		}
	}
	if (procID == Proc_Pelastic) {
		Er = f_chipXsection[ii]->GetRandom(Pthr + Pbinding + Mn, Tr_max + Mn);
	} else if (procID == Proc_Eelastic) {
		Er = f_chieXsection[ii]->GetRandom(Ethr + Me, Tr_max + Me);
	} else if (procID == Proc_Nuclelastic) {
		Er = f_chinuclXsection[ii]->GetRandom(Nuclthr, Tr_max); //Here the variable is the KINETIC energy of the recoiling nucleus
		Er = Er + Mnucl;
	}
	/*1a: correct the proton energy for binding effects*/
	if (procID == Proc_Pelastic) {
		Er = Er - Pbinding; /*Effective binding energy correction*/
	}

	/*1b: compute x-section total . No time consuming, since integrals are cached!*/
	if (procID == Proc_Pelastic) {
		sigma = f_chipXsection[ii]->Integral(Pthr + Pbinding + Mn, Tr_max + Mn);
	} else if (procID == Proc_Eelastic) {
		sigma = f_chieXsection[ii]->Integral(Ethr + Me, Tr_max + Me);
	} else if (procID == Proc_Nuclelastic) {
		sigma = f_chinuclXsection[ii]->Integral(Nuclthr, Tr_max); //Has to be integrated in this range since the variable is the KINETIC energy here (and not the total as before)
	}

	/*2: compute recoil chi TOTAL energy*/
	if (procID == Proc_Pelastic) {
		Echi = E0 + Mn - Er;
	} else if (procID == Proc_Eelastic) {
		Echi = E0 + Me - Er;
	} else if (procID == Proc_Nuclelastic) {
		Echi = E0 + Mnucl - Er;
	}

	/*3: compute the momenta*/
	Pchi = sqrt(Echi * Echi - Mchi * Mchi);
	if (procID == Proc_Pelastic) {
		Pr = sqrt(Er * Er - Mn * Mn);
	} else if (procID == Proc_Eelastic) {
		Pr = sqrt(Er * Er - Me * Me);
	} else if (procID == Proc_Nuclelastic) {
		Pr = sqrt(Er * Er - Mnucl * Mnucl);
	}
	/*4: compute the angle of the recoil nucleon wrt the initial chi momentum direction*/
	if (procID == Proc_Pelastic) {
		ctheta_r = E0 * E0 - Echi * Echi + Er * Er - Mn * Mn;
		ctheta_r /= 2 * P0 * Pr;
	} else if (procID == Proc_Eelastic) {
		ctheta_r = E0 * E0 - Echi * Echi + Er * Er - Me * Me;
		ctheta_r /= 2 * P0 * Pr;
	} else if (procID == Proc_Nuclelastic) {
		ctheta_r = E0 * E0 - Echi * Echi + Er * Er - Mnucl * Mnucl;
		ctheta_r /= 2 * P0 * Pr;
	}
	if (ctheta_r > 1) ctheta_r = 1;
	if (ctheta_r < -1) ctheta_r = -1;
	stheta_r = sqrt(1 - ctheta_r * ctheta_r);
	/*5: The azimuthal angle (around the incoming chi momentum direction) is flat*/
	phi_r = Rand.Uniform(-PI, PI);

	/*6: Now set the 4-vectors*/
	/*6a: build an orthogonal coordinate system, with v0 along the initial chi momentum direction*/
	v0 = chi.Vect().Unit();
	v1 = v0.Orthogonal();
	v1 = v1.Unit();
	v2 = v0.Cross(v1); //v2 = v0 x v1

	/*write the 3-momenta*/
	pr = v0 * Pr * ctheta_r + v1 * Pr * stheta_r * sin(phi_r) + v2 * Pr * stheta_r * cos(phi_r);
	pchi = p0 - pr;

	/*6b:Set them */
	recoil.SetVect(pr);
	recoil.SetE(Er);
	recoil_chi.SetVect(pchi);
	recoil_chi.SetE(Echi);

	return sigma;

}

/*Given the impact point on the front face (vin) and the incoming particle LorentzVector (chi for invisible decay, A' for visible),
 *  determine the interaction point within the fiducial volume and save it in vhit.
 Use a random distribution along the chi flight path, with uniform probability
 This function returns the length (in m) of the trajectory within the fiducial volume.
 Displacement is the lateral displacement (in m) of the detector, along x
 */
double KinUtils::findInteractionPoint(const TLorentzVector &chi, const TVector3 &fiducialV, const TVector3 &vin, TVector3 &vout, TVector3 &vhit) {

	double tz, tx, ty, tout, L;
	int sigPx, sigPy;

	sigPx = (chi.Px() > 0 ? 1 : -1);
	sigPy = (chi.Py() > 0 ? 1 : -1);

	tz = fiducialV.Z() / chi.Pz();
	tx = (sigPx * fiducialV.X() / 2 - vin.X()) / chi.Px();
	ty = (sigPy * fiducialV.Y() / 2 - vin.Y()) / chi.Py();
	tout = 0;

	if ((tz < tx) && (tz < ty)) {
		tout = tz;
	} else if ((tx < ty) && (tx < tz)) {
		tout = tx;
	} else if ((ty < tx) && (ty < tz)) {
		tout = ty;
	}
	vout.SetXYZ(tout * chi.Px() + vin.X(), tout * chi.Py() + vin.Y(), tout * chi.Pz() + vin.Z());
	vhit.SetXYZ(Rand.Uniform(vin.X(), vout.X()), Rand.Uniform(vin.Y(), vout.Y()), Rand.Uniform(vin.Z(), vout.Z()));
	L = (vout - vin).Mag();

	return L;
}

//This is the function to determine the interaction point in case of a cylinder (case1), for BDXmini.
//The cylinder is with the axis along y(vertical), center at x=0,y=0,z=ldet, radius is R,height is h
double KinUtils::findInteractionPointCylinder1(const TLorentzVector &chi,double ldet,double h,double R,TVector3 &vin,TVector3 &vout,TVector3 &vhit){

	double tIN,tOUT,t2,t3,t,L;

	double px=chi.Px();
	double py=chi.Py();
	double pz=chi.Pz();

	double delta=pz*pz*ldet*ldet-(pz*pz+px*px)*(ldet*ldet-R*R);
	if (delta<0){
		cout<<"KinUtils::findInteractionPointCylinder1 error, the delta is: "<<delta<<endl;
		return 0;
	}

	//entry point
	tIN = (pz*ldet - sqrt(delta))/(px*px+pz*pz);
	t2 = (pz*ldet + sqrt(delta))/(px*px+pz*pz);

	t3 = (h/2)/py;


	if ((t3>0)&&(t3<t2)&&(t3>tIN)){
		tOUT=t3;
	}else{
		tOUT=t2;
	}

	t=Rand.Uniform(tIN,tOUT);

	vin.SetXYZ(tIN*px,tIN*py,tIN*pz);
	vout.SetXYZ(tOUT*px,tOUT*py,tOUT*pz);

	vhit.SetXYZ(t*px,t*py,t*pz);


	L = (vout - vin).Mag();



	return L;
}

bool KinUtils::intersectsCylinder1(const TLorentzVector &chi,double ldet,double h,double R){

	double t1,y;

	double px=chi.Px();
	double py=chi.Py();
	double pz=chi.Pz();

	double delta=pz*pz*ldet*ldet-(pz*pz+px*px)*(ldet*ldet-R*R);
	if (delta<0) return false;

	t1 = (pz*ldet - sqrt(delta))/(px*px+pz*pz); //this is where it enters

	y = t1*py;
	if ((y>h/2)||(y<-h/2)) return false;

	return true;


}



