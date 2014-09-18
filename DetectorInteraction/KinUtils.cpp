#include <iostream>
#include <cmath>
#include "KinUtils.h"


using namespace std;

KinUtils::KinUtils(const double &m_Ebeam,const double &m_Mchi,const double &m_Maprime,const double &m_Msplit,const int &m_Seed):
Ebeam(m_Ebeam),
Mchi(m_Mchi),
Maprime(m_Maprime),
Msplit(m_Msplit),
Seed(m_Seed),
eTHR(0),
pTHR(0){

		Rand.SetSeed(m_Seed);
		f_chipXsection = new TF1*[100];
		f_chieXsection = new TF1*[100];

		for (int ii = 0; ii < 100; ii++) {
			f_chipXsection[ii] = new TF1(Form("f_chipXsection_%i", ii),this,&KinUtils::Er_chipXsection, 0, Ebeam, 2);
			f_chipXsection[ii]->SetNpx(1000);
			f_chipXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / 100);
			f_chieXsection[ii] = new TF1(Form("f_chieXsection_%i", ii),this,&KinUtils::Er_chieXsection, 0, Ebeam, 2);
			f_chieXsection[ii]->SetNpx(1000);
			f_chieXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / 100);
		}



}



/*double Tr_chipXsection
 Returns the chi p -> chi p DIFFERENTIAL cross section, dsigma/dEr, where Er is the recoil proton TOTAL energy
 Parameters: x[0]:   Recoil proton total energy
             par[0]: E0, total energy of incoming chi
*/
double KinUtils::Er_chipXsection(double *x,double *par){
	double Er=x[0];
	double E0=par[0];


	return 1;

}

/*double Tr_chieXsection
 Returns the chi e -> chi e DIFFERENTIAL cross section, dsigma/dTr, where Tr is the recoil electron kinetic energy
 Parameters: x[0]: Tr
             par[0]: T0, kinetic energy of incoming chi
*/
double KinUtils::Er_chieXsection(double *x,double *par){
	double Tr=x[0];
	double T0=par[0];
	double Tchi=T0-Tr;

	return 1;
}



int KinUtils::doProtonRecoil(const TLorentzVector &chi,TLorentzVector &recoil_p,TLorentzVector &recoil_chi){


	TVector3 v0,v1,v2;
	TVector3 p0,pr,pchi;
	double T0,Tr,Tchi;
	double P0,Pr,Pchi;
	double ctheta_r,stheta_r,phi_r;
	int ii;
	T0=chi.E()-Mchi;
	P0=chi.P();
	p0=chi.Vect();

	/*1: extract the recoil kinetic energy from the cross-section*/
	ii=0;
	ii=(int)rint((Ebeam/100)/T0);
	if (ii>=100) ii=100;
	Tr=f_chipXsection[ii]->GetRandom(pTHR,T0);

        /*2: compute recoil chi kinetic energy*/
	Tchi=T0-Tr;
	/*3: compute the momenta*/
	Pchi=sqrt((Tchi+Mchi)*(Tchi+Mchi)-Mchi*Mchi);
	Pr=sqrt((Tr+Mn)*(Tr+Mn)-Mn*Mn);
        /*4: compute the angle of the recoil nucleon wrt the initial chi momentum direction*/
	ctheta_r=(P0*P0+Pr*Pr-Pchi*Pchi)/(2*P0*Pr);
	if (ctheta_r>1) ctheta_r=1;
	if (ctheta_r<-1) ctheta_r=-1;
	stheta_r=sqrt(1-ctheta_r*ctheta_r);
	/*5: The azimuthal angle (around the incoming chi momentum direction) is flat*/
	phi_r=Rand.Uniform(-PI,PI);

	/*6: Now set the 4-vectors*/
	/*6a: build an orthogonal coordinate system, with v0 along the initial chi momentum direction*/
	v0=chi.Vect().Unit();
	v1=v0.Orthogonal();
	v1=v1.Unit();
	v2=v0.Cross(v1); //v2 = v0 x v1

	/*write the 3-momenta*/
	pr=v0*Pr*ctheta_r+v1*Pr*stheta_r*sin(phi_r)+v2*Pr*stheta_r*cos(phi_r);
	pchi=p0-pr;

	/*6b:Set them */
	recoil_p.SetVect(pr);
	recoil_p.SetE(Tr+Mn);
	recoil_chi.SetVect(pchi);
	recoil_chi.SetE(Tchi+Mchi);




	return 1;

}

int KinUtils::doElectronRecoil(const TLorentzVector &chi,TLorentzVector &recoil_e,TLorentzVector &recoil_chi){


	TVector3 v0,v1,v2;
	TVector3 p0,pr,pchi;
	double T0,Tr,Tchi;
	double P0,Pr,Pchi;
	double ctheta_r,stheta_r,phi_r;
	T0=chi.E()-Mchi;
	P0=chi.P();
	p0=chi.Vect();

	int ii;
	/*1: extract the recoil kinetic energy from the cross-section*/
	ii=0;
	ii=(int)rint((Ebeam/100)/T0);
	if (ii>=100) ii=100;
	Tr=f_chieXsection[ii]->GetRandom(eTHR,T0);
        /*2: compute recoil chi kinetic energy*/
	Tchi=T0-Tr;
	/*3: compute the momenta*/
	Pchi=sqrt((Tchi+Mchi)*(Tchi+Mchi)-Mchi*Mchi);
	Pr=sqrt((Tr+Me)*(Tr+Me)-Me*Me);
        /*4: compute the angle of the recoil nucleon wrt the initial chi momentum direction*/
	ctheta_r=(P0*P0+Pr*Pr-Pchi*Pchi)/(2*P0*Pr);
	if (ctheta_r>1) ctheta_r=1;
	if (ctheta_r<-1) ctheta_r=-1;
	stheta_r=sqrt(1-ctheta_r*ctheta_r);
	/*5: The azimuthal angle (around the incoming chi momentum direction) is flat*/
	phi_r=Rand.Uniform(-PI,PI);

	/*6: Now set the 4-vectors*/
	/*6a: build an orthogonal coordinate system, with v0 along the initial chi momentum direction*/
	v0=chi.Vect().Unit();
	v1=v0.Orthogonal();
	v1=v1.Unit();
	v2=v0.Cross(v1); //v2 = v0 x v1

	/*write the 3-momenta*/
	pr=v0*Pr*ctheta_r+v1*Pr*stheta_r*sin(phi_r)+v2*Pr*stheta_r*cos(phi_r);
	pchi=p0-pr;

	/*6b:Set them */
	recoil_e.SetVect(pr);
	recoil_e.SetE(Tr+Me);
	recoil_chi.SetVect(pchi);
	recoil_chi.SetE(Tchi+Mchi);

	return 1;
}













/*Given the impact point on the front face (vin) and the incoming particle LorentzVector (chi for invisible decay, A' for visible), determine the interaction point within the fiducial volume vhit.
Use a random distribution along the chi flight path, with uniform probability
*/
int KinUtils::findInteractionPoint(const TLorentzVector &chi,const TVector3 &fiducialV,const TVector3 &vin,TVector3 &vhit){

	TVector3 vout; //the exit point from the fiducial volume
	double tz,tx,ty,tout;
	int sigPx,sigPy;

	sigPx = ( chi.Px()>0 ? 1 : -1);
        sigPy = ( chi.Py()>0 ? 1 : -1);

	tz=fiducialV.Z()/chi.Pz();
	tx=(sigPx*fiducialV.X()/2-vin.X())/chi.Px();
	ty=(sigPy*fiducialV.Y()/2-vin.Y())/chi.Py();
	tout=0;
	if ((tz<tx) && (tz<ty)){
		tout=tz;
	}
	else if ((tx<ty)&&(tx<tz)){
		tout=tx;
	}
	else if ((ty<tx)&&(ty<tz)){
		tout=ty;
	}
	vout.SetXYZ(tout*chi.Px()+vin.X(),tout*chi.Py()+vin.Y(),tout*chi.Pz());
	vhit.SetXYZ(Rand.Uniform(vin.X(),vout.X()),Rand.Uniform(vin.Y(),vout.Y()),Rand.Uniform(vin.Z(),vin.Z()+vout.Z()));//along z shift wrt beam-dump
	return 1;
}
