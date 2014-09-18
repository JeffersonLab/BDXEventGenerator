#ifndef KinUtils_G
#define KinUtils_G



#define PI 3.14156
#define Me 0.511E-3
#define Mn 0.9383

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TRandom3.h"


class KinUtils{

	private:
		TF1** f_chipXsection;
		TF1** f_chieXsection;

		double Ebeam;
		double Mchi,Maprime,Msplit; //chi and aprime mass in GeV
		double pTHR,eTHR; //the two thresholds in GeV for the recoil proton kinetic energy and the recoil electron kinetic energy.
		int Seed;
		TRandom3 Rand;

	public:
		KinUtils(const double &m_Ebeam,const double &m_Mchi,const double &m_Maprime,const double &m_Msplit,const int &m_Seed);
		int doProtonRecoil(const TLorentzVector &chi,TLorentzVector &proton,TLorentzVector &chiPrime);
		int doElectronRecoil(const TLorentzVector &chi,TLorentzVector &electron,TLorentzVector &chiPrime);
		int findInteractionPoint(const TLorentzVector &chi,const TVector3 &fiducialV,const TVector3 &vin,TVector3 &vhit);
		double Er_chipXsection(double *x,double *par);
		double Er_chieXsection(double *x,double *par);


};






#endif
