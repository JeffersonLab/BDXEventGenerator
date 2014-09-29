#include <iostream>
#include <cmath>
#include "KinUtils.h"


using namespace std;

KinUtils::KinUtils(const double &m_Ebeam,const double &m_Mchi,const double &m_Maprime,const double &m_Msplit,const double &m_Epsilon,const double &m_AlphaD,const int &m_Seed):
Ebeam(m_Ebeam),
Mchi(m_Mchi),
Maprime(m_Maprime),
Msplit(m_Msplit),
Seed(m_Seed),
Ethr(0),
Pthr(0),
Epsilon(m_Epsilon),
AlphaDark(m_AlphaD),
Alpha(1./137.){

		Rand.SetSeed(m_Seed);
		f_chipXsection = new TF1*[100];
		f_chieXsection = new TF1*[100];

		/*This part is really critical. The cross-section is a function of the final state recoil energy, and as "parameter" needs the incoming
		 * chi energy. However, if event-by-event we set this parameter, and then call a "GetRandom", we trigger the integration computation for ALL events,
		 * and this is very time-consuming.
		 * Instead, I "bin" wrt the incoming chi energy (100 bins), and then extract the random number from the function defining the energy in that bin!
		 */
		cout<<"KinUtils::KinUtils setting the cross-section functions";
		for (int ii = 0; ii < 100; ii++) {
			f_chipXsection[ii] = new TF1(Form("f_chipXsection_%i", ii),this,&KinUtils::Er_chipXsection, Pthr+Mn, Ebeam+Mn-Mchi, 1);
			f_chipXsection[ii]->SetNpx(1000);
			f_chipXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / 100);
			f_chieXsection[ii] = new TF1(Form("f_chieXsection_%i", ii),this,&KinUtils::Er_chieXsection, Ethr+Me, Ebeam+Me-Mchi, 1);
			f_chieXsection[ii]->SetNpx(1000);
			f_chieXsection[ii]->FixParameter(0, (ii + 1) * Ebeam / 100);
			//cout<<f_chipXsection[ii]->Integral(Pthr+Mn,Ebeam)<<" "<<f_chieXsection[ii]->Integral(Ethr+Me,Ebeam)<<endl;
			f_chipXsection[ii]->Integral(Pthr+Mn,Ebeam+Mn-Mchi);
			f_chieXsection[ii]->Integral(Ethr+Me,Ebeam+Mn-Mchi);
		}

		cout<<"KinUtils::KinUtils created"<<endl;
		PrintParameters();


}

void KinUtils::PrintParameters(){
			cout<<"KinUtils Parameters:"<<endl;
			cout<<"e- beam (GeV): \t\t"<<Ebeam<<endl;
			cout<<"Maprime (GeV): \t \t"<<Maprime<<endl;
			cout<<"Mchi (GeV): \t \t"<<Mchi<<endl;
			cout<<"Split (GeV): \t  \t"<<Msplit<<endl;
			cout<<"Epsilon: \t \t"<<Epsilon<<endl;
			cout<<"AlphaDark: \t \t"<<AlphaDark<<endl;
			cout<<"e- threshold: \t \t"<<Ethr<<endl;
			cout<<"p threshold: \t \t"<<Pthr<<endl;
}

/*double Tr_chipXsection
 Returns the chi p -> chi p DIFFERENTIAL cross section, dsigma/dEr, where Er is the recoil proton TOTAL energy
 Parameters: x[0]:   Recoil proton total energy
             par[0]: E0, total energy of incoming chi
*/
double KinUtils::Er_chipXsection(double *x,double *par){
	double Er=x[0];
	double E0=par[0];


	//The reaction is: chi(p1)+N(p2)->chi(k1)+N(k2).
	//Perform the cross-section in the lab frame, where p1=(0,0,P,E), p2=(0,0,0,Mn), k2=(xx,yy,zz,Er)
	//1: Compute Lorentz-invariant dot products

	double p1p1=Mchi*Mchi;
	double k1k1=Mchi*Mchi;
	double p1p2=E0*Mn;
	double p1k1=-0.5*(Mn*Mn-Mchi*Mchi+Er*Mchi); //A.C. Ask Eder to check
	double p1k2=Mn*(E0+Mn-Er);
	double p2k1=(E0+Mn-Er)*Mn;
	double p2k2=Er*Mn;
	double k1k2=Mn*Er;

	double k1p2=p2k1;
	double k2p2=p2k2;
	double k1p1=p1k1;

	double t=2*Mn*Mn-2*Mn*Er;
	double s=Mn*Mn + Mchi*Mchi + 2*Mn*E0;

	//2: the amplitude squared
	double  ampsq=1.0*(k1k2*p1p2 + p1k2*k1p2 - Mchi*Mchi*k2p2 - Mn*Mn*k1p1 + 2*Mchi*Mchi*Mn*Mn)/(pow((t-Maprime*Maprime),2));
	if (ampsq<0.0){
	            cout<<"Er_chipXsection: negative amplitude!!!"<<endl;
	            return 0;
	}
	//3: the 2 momenta in the CM frame (before/after). Since this is elastic scattering, they're the same!!!

	double p=((pow(s-Mchi*Mchi - Mn*Mn,2) - 4*Mchi*Mchi*Mn*Mn)/(4*s));
	p=sqrt(p);
	double k=p;

	//4: the Jacobian for the transformation cos(theta_CM) --> Recoil total energy in LAB frame
	double dcostheta = Mn / (p*k); //AC: checked this formula.

	//5: the phase-space
	double S=(k/s*p)*dcostheta;

	 //7: the Form-Factor, using a dipole approximation.
	double FF=1/pow(1+(Er*Er-Mn*Mn)/Mn*Mn,2);

	//8: put everything together

	double dsigma=pow(FF,2)*ampsq*S; //in "GeV^-2 units" (and no coupling yet)
	dsigma = dsigma * Epsilon*Epsilon*4*PI*Alpha*AlphaDark;
	dsigma *=GeVm2cm2;



	return dsigma;

}

/*double Tr_chieXsection
 Returns the chi e -> chi e DIFFERENTIAL cross section, dsigma/dTr, where Tr is the recoil electron kinetic energy
 Parameters: x[0]: Tr
             par[0]: T0, kinetic energy of incoming chi
*/
double KinUtils::Er_chieXsection(double *x,double *par){
	double Er=x[0];
	double E0=par[0];

	//The reaction is: chi(p1)+e(p2)->chi(k1)+e(k2).
	//Perform the cross-section in the lab frame, where p1=(0,0,P,E), p2=(0,0,0,Me), k2=(xx,yy,zz,Er)
	//1: Compute Lorentz-invariant dot products

	double p1p1=Mchi*Mchi;
	double k1k1=Mchi*Mchi;
	double p1p2=E0*Me;
	double p1k1=-0.5*(Me*Me-Mchi*Mchi+Er*Mchi);
	double p1k2=Me*(E0+Me-Er);
	double p2k1=(E0+Me-Er)*Me;
	double p2k2=Er*Me;
	double k1k2=Me*Er;

	double k1p2=p2k1;
	double k2p2=p2k2;
	double k1p1=p1k1;

	double t=2*Me*Me-2*Me*Er;
	double s=Me*Me + Mchi*Mchi + 2*Me*E0;

	//2: the amplitude squared
	double  ampsq=1.0*(k1k2*p1p2 + p1k2*k1p2 - Mchi*Mchi*k2p2 - Me*Me*k1p1 + 2*Mchi*Mchi*Me*Me)/(pow((t-Maprime*Maprime),2));
	if (ampsq<0.0){
	            cout<<"Er_chieXsection: negative amplitude!!! "<<E0<<" "<<Er<<endl;
	            return 0;
	}
	//3: the 2 momenta in the CM frame (before/after). Since this is elastic scattering, they're the same!!!
	double p=((pow(s-Mchi*Mchi - Me*Me,2) - 4*Mchi*Mchi*Me*Me)/(4*s));
	p=sqrt(p);
	double k=p;
	//4: the Jacobian for the transformation cos(theta_CM) --> Recoil total energy in LAB frame
	double dcostheta = Me / (p*k); //AC: checked this formula.

	//5: the phase-space
	double S=(k/s*p)*dcostheta;

	//8: put everything together
	double dsigma=ampsq*S; //in "GeV^-2 units" (and no coupling yet)
	dsigma = dsigma * Epsilon*Epsilon*4*PI*Alpha*AlphaDark;
	dsigma *=GeVm2cm2;
	return dsigma;
}


/*This is the routine that, given the cross-section for the ELASTIC interaction chi-e --> chi-e OR chi-p --> chi-p
 * estracts the recoil kinematics.
 * procID is an integer, for the process under study:
 *
 * Returns: the total cross-section, in cm2, for the interaction process, integrated over final state (with the kinetic energy cut)
 */

double KinUtils::doElasticRecoil(const TLorentzVector &chi,TLorentzVector &recoil,TLorentzVector &recoil_chi,const int &procID){


	TVector3 v0,v1,v2;
	TVector3 p0,pr,pchi;
	double E0,Er,Echi;
	double P0,Pr,Pchi;
	double ctheta_r,stheta_r,phi_r,ret;
	int ii;
	E0=chi.E();
	P0=chi.P();
	p0=chi.Vect();

	/*1: extract the recoil total energy from the cross-section*/
	ii=0;
	ii=(int)(E0/(Ebeam/100));
	if (ii>=100) ii=99; //should not happen

	(procID==Proc_Pelastic ? Er=f_chipXsection[ii]->GetRandom(Pthr+Mn,E0+Mn-Mchi):Er=f_chieXsection[ii]->GetRandom(Ethr+Mn,E0+Me-Mchi));

	/*2: compute recoil chi total energy*/
	Echi=E0+Mn-Er;
	/*3: compute the momenta*/
	Pchi=sqrt(Echi*Echi-Mchi*Mchi);
	(procID==Proc_Pelastic ? Pr=sqrt(Er*Er-Mn*Mn) : Pr=sqrt(Er*Er-Me*Me));
    /*4: compute the angle of the recoil nucleon wrt the initial chi momentum direction*/
	ctheta_r=(P0*P0+Pr*Pr-Pchi*Pchi)/(2*P0*Pr);
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
	recoil.SetVect(pr);
	recoil.SetE(Er);
	recoil_chi.SetVect(pchi);
	recoil_chi.SetE(Echi);

	//no time consuming, since integrals are cached!
	(procID==Proc_Pelastic ? ret=f_chipXsection[ii]->Integral(Pthr+Mn,E0+Mn-Mchi):ret=f_chieXsection[ii]->Integral(Ethr+Mn,E0+Me-Mchi));
	return ret;


}












/*Given the impact point on the front face (vin) and the incoming particle LorentzVector (chi for invisible decay, A' for visible),
 *  determine the interaction point within the fiducial volume and save it in vhit.
Use a random distribution along the chi flight path, with uniform probability
This function returns the length (in m) of the trajectory within the fiducial volume.
*/
double KinUtils::findInteractionPoint(const TLorentzVector &chi,const TVector3 &fiducialV,const TVector3 &vin,TVector3 &vhit){

	TVector3 vout; //the exit point from the fiducial volume
	double tz,tx,ty,tout,L;
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
	vout.SetXYZ(tout*chi.Px()+vin.X(),tout*chi.Py()+vin.Y(),tout*chi.Pz()+vin.Z());
	vhit.SetXYZ(Rand.Uniform(vin.X(),vout.X()),Rand.Uniform(vin.Y(),vout.Y()),Rand.Uniform(vin.Z(),vin.Z()+vout.Z()));//along z shift wrt beam-dump;
	L=(vout-vin).Mag();




	return L;
}
