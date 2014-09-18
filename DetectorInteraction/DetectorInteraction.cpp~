
#include <iostream>
#include <map>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TChain.h"
#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"

#include "LHEF.h"

#include "ExRootAnalysis/ExRootClasses.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#define Me 0.511E-3
#define Mn 0.9383
#define PI 3.1415


using namespace std;


TF1 **f_chipXsection;
TF1 **f_chieXsection;

int doProtonRecoil(const TLorentzVector &chi,TLorentzVector &proton,TLorentzVector &chiPrime);
int doElectronRecoil(const TLorentzVector &chi,TLorentzVector &electron,TLorentzVector &chiPrime);
int findInteractionPoint(const TLorentzVector &chi,const TVector3 &fiducialV,const TVector3 &vin,TVector3 &vhit);
double Tr_chipXsection(double *x,double *par);
double Tr_chieXsection(double *x,double *par);

double Ebeam; //the primary beam energy, i.e. the maximum chi energy in lab frame GeV
double m_chi,m_aprime; //chi and aprime mass in GeV
double pTHR,eTHR; //the two thresholds in GeV for the recoil proton kinetic energy and the recoil electron kinetic energy.

TRandom3 m_rand;
//---------------------------------------------------------------------------

void AnalyseParticles(LHEF::Reader *reader)
{
  LHEF::HEPEUP &hepeup = reader->hepeup;
  const LHEF::HEPRUP &heprup = reader->heprup;
  
  Int_t particle, n_inside;
  long PID;
  Double_t signPz, cosTheta, M;
 
  Double_t xthetaxfmin,xthetaxfmax,xthetayfmin,xthetayfmax;
  TLorentzVector chi,recoil_chi,recoil_p,recoil_e;
  TVector3 vin,vhit,fiducialV;
  vector<double> tmp_v;	
	
	
  n_inside=0;
  xthetaxfmin=xthetayfmin=0;
  
  xthetaxfmax=(heprup.lx/(2*heprup.ldet)); //ok without ATAN for this check
  xthetayfmax=(heprup.ly/(2*heprup.ldet)); //ok without ATAN for this check
  
  fiducialV.SetXYZ(heprup.lx,heprup.ly,heprup.lz);
	
 for(particle = 0; particle < hepeup.NUP; ++particle)
  {
	  
    PID = hepeup.IDUP[particle];    
    if ((PID!=-611)&&(PID!=611)) continue; //other particles are ok. Go on
    
    M = hepeup.PUP[particle][4];
    chi.SetPxPyPzE(hepeup.PUP[particle][0],hepeup.PUP[particle][1],hepeup.PUP[particle][2],hepeup.PUP[particle][3]);
    cosTheta = TMath::Abs(chi.CosTheta());
    signPz = (chi.Pz() >= 0.0) ? 1.0 : -1.0;
    /*In MADGRAPH, the fiducial volume cut was inclusive, i.e. the event was considered "good" if at least one of the two chi were inside.
	Now, I need to find which of the two are inside, and use that for the interaction.
	If both are inside, I take only one.
	If none is inside, there is an error!*/
    if ((fabs(chi.Px()/chi.Pz()) < xthetaxfmax) && (fabs(chi.Py()/chi.Pz()) < xthetayfmax)) n_inside++;	  
    else{ //for now, just mark this chi out of the fiducial volume with a status "0"
	    hepeup.ISTUP[particle]=0;
	    continue;
    }
    //now check if this is the second chi, if it is within the fiducial volume, and if the first was already in the fiducial volume. If so, we skip the second chi 
    if (n_inside==2){
	    hepeup.ISTUP[particle]=0;  
	    continue;
    }
   vin.SetXYZ((chi.Px()/chi.Pz())*heprup.ldet,(chi.Py()/chi.Pz())*heprup.ldet,heprup.ldet); //the chi hit position in the fiducial volume front-face	
    //From here, we have a chi within the fiducial volume.	
   // See which interaction to consider
  	  
   switch (heprup.procid){
   case 0: //nothing to do
	   break;
   case 1: //proton elastic	   
	   doProtonRecoil(chi,recoil_p,recoil_chi);
	   findInteractionPoint(chi,fiducialV,vin,vhit);
	   //add particles to hepeup
	   //final state chi
	   hepeup.IDUP.push_back(9611);
	   hepeup.ISTUP.push_back(1);
	   hepeup.MOTHUP.push_back(std::make_pair(particle+1,particle+1));
           hepeup.ICOLUP.push_back(std::make_pair(0,0));
           tmp_v.clear();
	   tmp_v.push_back(recoil_chi.Px());
	   tmp_v.push_back(recoil_chi.Py());
	   tmp_v.push_back(recoil_chi.Pz());
	   tmp_v.push_back(recoil_chi.E());
	   tmp_v.push_back(recoil_chi.M());	                          
           hepeup.PUP.push_back(tmp_v);
           hepeup.VTIMUP.push_back(0);
	   hepeup.SPINUP.push_back(0);
           //final state proton                        
	   hepeup.IDUP.push_back(92212);
	   hepeup.ISTUP.push_back(1);
	   hepeup.MOTHUP.push_back(std::make_pair(particle+1,particle+1));
	   hepeup.ICOLUP.push_back(std::make_pair(0,0));
	   tmp_v.clear();
	   tmp_v.push_back(recoil_p.Px());
	   tmp_v.push_back(recoil_p.Py());
	   tmp_v.push_back(recoil_p.Pz());
	   tmp_v.push_back(recoil_p.E());
	   tmp_v.push_back(recoil_p.M());	                          
	   hepeup.PUP.push_back(tmp_v);
	   hepeup.VTIMUP.push_back(0);
	   hepeup.SPINUP.push_back(0);       
	   hepeup.NUP+=2;
	   break;
   case 2: //electron elastic	  
	   doElectronRecoil(chi,recoil_e,recoil_chi);
	   findInteractionPoint(chi,fiducialV,vin,vhit);
	   //add particles to hepeup
	   //final state chi
	   hepeup.IDUP.push_back(9611);
	   hepeup.ISTUP.push_back(1);
	   hepeup.MOTHUP.push_back(std::make_pair(particle+1,particle+1));
	   hepeup.ICOLUP.push_back(std::make_pair(0,0));
	   tmp_v.clear();
	   tmp_v.push_back(recoil_chi.Px());
           tmp_v.push_back(recoil_chi.Py());
           tmp_v.push_back(recoil_chi.Pz());
           tmp_v.push_back(recoil_chi.E());
           tmp_v.push_back(recoil_chi.M());	                          
           hepeup.PUP.push_back(tmp_v);
           hepeup.VTIMUP.push_back(0);
           hepeup.SPINUP.push_back(0);
           //final state electron                        
	   hepeup.IDUP.push_back(911);
           hepeup.ISTUP.push_back(1);
	   hepeup.MOTHUP.push_back(std::make_pair(particle+1,particle+1));
	   hepeup.ICOLUP.push_back(std::make_pair(0,0));
	   tmp_v.clear();
	   tmp_v.push_back(recoil_e.Px());
	   tmp_v.push_back(recoil_e.Py());
	   tmp_v.push_back(recoil_e.Pz());
	   tmp_v.push_back(recoil_e.E());
	   tmp_v.push_back(recoil_e.M());	                          
	   hepeup.PUP.push_back(tmp_v);
	   hepeup.VTIMUP.push_back(0);
	   hepeup.SPINUP.push_back(0);          
	   hepeup.NUP+=2;
	   break;
   case 3: //inelastic
	   
	   break;
   default:
	   cout<<"Error, interaction not recognized"<<endl;
	   break;
    }	
	  
	  //use the eventComments for the vertex location (in m, in the form x y z)
	  reader->eventComments=Form("%f %f %f",vhit.X(),vhit.Y(),vhit.Z());
  }	
	
  if (n_inside==0) cout<<"Error with this event, no chi inside fiducial volume (but should be at least one!)"<<endl;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "ExRootLHEFConverter";

  if(argc != 3)
  {
    cout << " Usage: " << appName << " input_file" << " output_file" << endl;
    cout << " input_file - input file in LHEF format," << endl;
    cout << " output_file - output file in LHEF format." << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  
	
  TApplication app(appName, &appargc, appargv);

  // Open a stream connected to an event file:
  ifstream inputFileStream(argv[1]);
  ofstream outputFileStream(argv[2]);
	
  // Create the Reader object:
	LHEF::Reader *inputReader = new LHEF::Reader(inputFileStream); //this triggers also "init"
  // Create the Writer object:
  LHEF::Writer *outputWriter = new LHEF::Writer(outputFileStream);
  outputWriter->heprup=inputReader->heprup;
  outputWriter->headerStream.str(inputReader->headerBlock);
  outputWriter->init();	
	
	//instear of having, for each process, a single function with a free parameter (the incident chi energy),
	//it is much better to have 100 of them, each for a certain energy: otherwise, the integral gets calculated for each event!
        cout<< " ** Preparing the cross-section functions ** " <<endl;
	Ebeam=inputReader->heprup.EBEAM;
	m_chi=inputReader->heprup.FMASS;
	m_aprime=inputReader->heprup.APMASS;
	m_rand.SetSeed(inputReader->heprup.SEED);
	f_chipXsection=new TF1*[100];
        f_chieXsection=new TF1*[100];
	
	
   for (int ii=0;ii<100;ii++){
 	f_chipXsection[ii]=new TF1(Form("f_chipXsection_%i",ii),Tr_chipXsection,0,Ebeam,2);
 	f_chipXsection[ii]->SetNpx(1000);
	f_chipXsection[ii]->FixParameter(0,(ii+1)*Ebeam/100);
		
	f_chieXsection[ii]=new TF1(Form("f_chieXsection_%i",ii),Tr_chieXsection,0,Ebeam,2);
	f_chieXsection[ii]->SetNpx(1000);
	f_chieXsection[ii]->FixParameter(0,(ii+1)*Ebeam/100);
   }
	
	
  cout << "** Calculating number of events to process. Please wait..." << endl;
  Long64_t allEntries = inputReader->getNumberOfEvents();
  cout << "** Input file contains " << allEntries << " events" << endl;

  if(allEntries > 0)
  {
    ExRootProgressBar progressBar(allEntries);
    
    // Loop over all events
    Long64_t entry = 0;
    while(inputReader->readEvent())
    {
     
      //This is the function that triggers the interaction in the fiducial volume.
      if (inputReader->heprup.procid) AnalyseParticles(inputReader);
      
      outputWriter->hepeup=inputReader->hepeup;
      outputWriter->eventStream.str(inputReader->eventComments);
      outputWriter->writeEvent();
    

      progressBar.Update(entry);

      ++entry;
    }

    progressBar.Finish();
  }


  cout << "** Exiting..." << endl;

  
  delete inputReader;
  delete outputWriter;
}









/*double Tr_chipXsection
 Returns the chi p -> chi p DIFFERENTIAL cross section, dsigma/dTr, where Tr is the recoil proton kinetic energy
 Parameters: x[0]:   Recoil proton kinetic energy
             par[0]: T0, kinetic energy of incoming chi
*/
double Tr_chipXsection(double *x,double *par){
	double Tr=x[0];
	double T0=par[0];
	double Tchi=T0-Tr;
	double NUM,DEN;
	double A,B,F1,F2;
	
	F1=1./((1+2*Tr/Mn)*(1+2*Tr/Mn));
	F2=1.79/((1+2*Tr/Mn)*(1+2*Tr/Mn));
	
	/*
	A=2*Mn*T0*(T0-Tr)-my_vars.m_chi*my_vars.m_chi*Tr;
	B=-Tr*((2*T0-Tr)*(2*T0-Tr)+2*Mn*Tr-4*my_vars.m_chi*my_vars.m_chi);
	F1=1./((1+2*Tr/Mn)*(1+2*Tr/Mn));
	F2=1.79/((1+2*Tr/Mn)*(1+2*Tr/Mn));
	NUM=4*PI*ALPHA*my_vars.alpha_dark*my_vars.epsilon_aprime*my_vars.epsilon_aprime;
	NUM*=(A*F1*F1-B*F2*F2*0.25);
	DEN=(my_vars.m_aprime*my_vars.m_aprime+2*Mn*Tr)*(my_vars.m_aprime*my_vars.m_aprime+2*Mn*Tr)*(T0*T0-my_vars.m_chi*my_vars.m_chi);	
	return N*NUM/DEN;
	*/
	
	NUM=F1*F1+(Tr/(2*Mn))*F2*F2;
	
	DEN=(m_aprime*m_aprime+2*Mn*Tr)*(m_aprime*m_aprime+2*Mn*Tr);
	
	return (2*Mn)*(NUM/DEN);
	
	
	
}

/*double Tr_chieXsection
 Returns the chi e -> chi e DIFFERENTIAL cross section, dsigma/dTr, where Tr is the recoil electron kinetic energy
 Parameters: x[0]: Tr
             par[0]: T0, kinetic energy of incoming chi
*/
double Tr_chieXsection(double *x,double *par){
	double Tr=x[0];
	double T0=par[0];
	double Tchi=T0-Tr;
	double NUM,DEN;
	

	NUM=(4*Me*m_chi*m_chi*Tr+pow(m_chi*m_chi+Me*(T0-Tr),2.));
	DEN=pow(m_aprime*m_aprime+2*Me*Tr,2.)*pow(m_chi*m_chi+2*Me*T0,2.);
	
	return (NUM/DEN);
}



int doProtonRecoil(const TLorentzVector &chi,TLorentzVector &recoil_p,TLorentzVector &recoil_chi){
	
	
	TVector3 v0,v1,v2;
	TVector3 p0,pr,pchi;
	double T0,Tr,Tchi;
	double P0,Pr,Pchi;		
	double ctheta_r,stheta_r,phi_r;
	int ii;
	T0=chi.E()-m_chi;
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
	Pchi=sqrt((Tchi+m_chi)*(Tchi+m_chi)-m_chi*m_chi);
	Pr=sqrt((Tr+Mn)*(Tr+Mn)-Mn*Mn);				
        /*4: compute the angle of the recoil nucleon wrt the initial chi momentum direction*/
	ctheta_r=(P0*P0+Pr*Pr-Pchi*Pchi)/(2*P0*Pr);
	if (ctheta_r>1) ctheta_r=1;
	if (ctheta_r<-1) ctheta_r=-1;
	stheta_r=sqrt(1-ctheta_r*ctheta_r);
	/*5: The azimuthal angle (around the incoming chi momentum direction) is flat*/
	phi_r=m_rand.Uniform(-PI,PI);	
	
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
	recoil_chi.SetE(Tchi+m_chi);
	
	
	
	
	return 1;
	
}

int doElectronRecoil(const TLorentzVector &chi,TLorentzVector &recoil_e,TLorentzVector &recoil_chi){
	
	
	TVector3 v0,v1,v2;
	TVector3 p0,pr,pchi;
	double T0,Tr,Tchi;
	double P0,Pr,Pchi;		
	double ctheta_r,stheta_r,phi_r;
	T0=chi.E()-m_chi;
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
	Pchi=sqrt((Tchi+m_chi)*(Tchi+m_chi)-m_chi*m_chi);
	Pr=sqrt((Tr+Me)*(Tr+Me)-Me*Me);				
        /*4: compute the angle of the recoil nucleon wrt the initial chi momentum direction*/
	ctheta_r=(P0*P0+Pr*Pr-Pchi*Pchi)/(2*P0*Pr);
	if (ctheta_r>1) ctheta_r=1;
	if (ctheta_r<-1) ctheta_r=-1;
	stheta_r=sqrt(1-ctheta_r*ctheta_r);
	/*5: The azimuthal angle (around the incoming chi momentum direction) is flat*/
	phi_r=m_rand.Uniform(-PI,PI);	
	
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
	recoil_chi.SetE(Tchi+m_chi);
	
	return 1;	
}













/*Given the impact point on the front face (vin) and the incoming particle LorentzVector (chi for invisible decay, A' for visible), determine the interaction point within the fiducial volume vhit.
Use a random distribution along the chi flight path, with uniform probability
*/
int findInteractionPoint(const TLorentzVector &chi,const TVector3 &fiducialV,const TVector3 &vin,TVector3 &vhit){
	
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
	vhit.SetXYZ(m_rand.Uniform(vin.X(),vout.X()),m_rand.Uniform(vin.Y(),vout.Y()),m_rand.Uniform(vin.Z(),vin.Z()+vout.Z()));//along z shift wrt beam-dump
	return 1;			
}









