#ifndef KinUtils_G
#define KinUtils_G



#define PI 3.14156
#define Me 0.511E-3
#define Mn 0.9383
#define GeVm2cm2 3.9204E-28       //This is the conversion factor for the cross section, from GeV^-2 to cm2

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TRandom3.h"

enum ProcID{
	Proc_nothing=0,
	Proc_Pelastic=1,
	Proc_Eelastic=2,
	Proc_Pinelastic=3,
	Proc_Einelastic=4
};

class KinUtils{

	private:
		TF1** f_chipXsection;
		TF1** f_chieXsection;

		double Ebeam;
		double Mchi,Maprime,Msplit; //chi and aprime mass in GeV
		double Epsilon,AlphaDark,Alpha;

		double Pthr,Ethr; //the two thresholds in GeV for the recoil proton kinetic energy and the recoil electron kinetic energy.
		int Seed;
		TRandom3 Rand;

		static const int nFunctionsElastic = 1000;


	public:


		KinUtils(const double &m_Ebeam,const double &m_Mchi,const double &m_Maprime,const double &m_Msplit,const double &m_Epsilon,const double &m_alphaD,const double &m_Pthr,const double &m_Ethr,const int &m_Seed);
		void PrintParameters();
		void setAlpha(double m_Alpha){Alpha=m_Alpha;};
		void setAlphaD(double m_AlphaDark){AlphaDark=m_AlphaDark;};
		void setEpsilon(double m_Epsilon){Epsilon=m_Epsilon;};

		double doElasticRecoil(const TLorentzVector &chi,TLorentzVector &recoil,TLorentzVector &chiPrime,const int &procID);
		double findInteractionPoint(const TLorentzVector &chi,const TVector3 &fiducialV,const TVector3 &vin,TVector3 &vout,TVector3 &vhit);
		double Er_chipXsection(double *x,double *par);
		double Er_chieXsection(double *x,double *par);







};






#endif
