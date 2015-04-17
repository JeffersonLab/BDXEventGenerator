#include <iostream>
#include "TClonesArray.h"
#include "TSystem.h"
#include "TChain.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

TRandom3 mrand(0);


//Acceptance at 2 m
double Nacc2a,Nacc2b,Nacc2c;  //20x20 cm2, 50x50 cm2, 100x100 cm2

//Acceptance at 5 m
double Nacc5a,Nacc5b,Nacc5c;  //20x20 cm2, 50x50 cm2, 100x100 cm2
	
//Acceptance at 10 m
double Nacc10a,Nacc10b,Nacc10c;    //20x20 cm2, 50x50 cm2, 100x100 cm2
	
//Acceptance at 15 m
double Nacc15a,Nacc15b,Nacc15c;   //20x20 cm2, 50x50 cm2, 100x100 cm2
	
//Acceptance at 20 m
double Nacc20a,Nacc20b,Nacc20c;   //20x20 cm2, 50x50 cm2, 100x100 cm2

void ana(string fname){
	
	
	TChain *chain = new TChain("LHEF");
	chain->Add(fname.c_str());
  //  chain->SetProof();
	int N=chain->GetEntries();
	cout<<N<<endl;
	
	
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchP = treeReader->UseBranch("Particle");
	TClonesArray *branchE = treeReader->UseBranch("Event");
	
	
	TRootLHEFEvent* evt;
	TRootLHEFParticle *particle,*chi,*achi;
	
	TLorentzVector Pchi,Pachi;
	TLorentzVector *Pchi_in;
	
	int Np;

	
	double x,y;
	double L,W;
	TVector3 vin,vout,vhit;
	
  //Histograms
	TH1D *hEnergy0=new TH1D("hEnergy0","hEnergy0; #chi Energy (GeV)",100,0,1.3);
	TH1D *hEnergy1=new TH1D("hEnergy1","hEnergy1",100,0,1.3);
	
	TH1D *hAngle0=new TH1D("hAngle0","hAngle0; #chi Angle (deg)",100,0,30);
	TH1D *hAngle1=new TH1D("hAngle1","hAngle1",100,0,30);
	
	TH2D *hSpace0=new TH2D("hSpace0","hSpace0; #chi X (at 1 m); #chi hit Y (at 1 m)",200,-2.,2.,200,-2.,2.);
	TH2D *hSpace1=new TH2D("hSpace1","hSpace1; #chi X (at 10 m); #chi hit Y (at 10 m)",200,-20.,20.,200,-20.,20.);
	
	TH2D *hVin0=new TH2D("hVin0","hVin0; #chi in X (m); #chi in Y (m)",100,-1.,1.,100,-1.,1.);
	TH2D *hVout0=new TH2D("hVout0","hVout0; #chi out X (m); #chi out Y (m)",100,-1.,1.,100,-1.,1.);
	
	TH1D *hL0=new TH1D("hL","hL; #chi length in det (m)",100,0.,5.); //all L
	TH1D *hL1=new TH1D("hL1","hL1",100,0.,5.); //L for particles exiting back
	TH1D *hL2=new TH1D("hL2","hL2",100,0.,5.); //L for particles exiting side
	
	TH2D *hLvsE0=new TH2D("hLvsE","hLvsE; E(GeV); L(m)",100,0.,1.3,100,0,5.);
	
	TH2D *hLvsTheta0=new TH2D("hLvsTheta0","hLvsTheta0; Theta(deg); L(m)",100,0.,20.,100,0,5.);
	
	Nacc2a=Nacc2b=Nacc2c=Nacc5a=Nacc5b=Nacc5c=Nacc10a=Nacc10b=Nacc10c=Nacc15a=Nacc15b=Nacc15c=Nacc20a=Nacc20b=Nacc20c=0;
	
	for (int ii=0;ii<N;ii++){
		treeReader->ReadEntry(ii);
		evt=(TRootLHEFEvent*)(branchE->At(0));
		
		Np=evt->Nparticles;
		L=evt->L;
		W=evt->W;
		vin.SetXYZ(evt->Vin.X(),evt->Vin.Y(),evt->Vin.Z());
		vout.SetXYZ(evt->Vout.X(),evt->Vout.Y(),evt->Vout.Z());
		vhit.SetXYZ(evt->Vhit.X(),evt->Vhit.Y(),evt->Vhit.Z());
		
		
		for (int jj=0;jj<Np;jj++){
			particle=(TRootLHEFParticle*)(branchP->At(jj));
			switch (particle->PID){
			case 611:
				chi=particle;
				break;
			case -611:
				achi=particle;
				break;
			default:
				break;
			}
		}
		Pchi.SetXYZT(chi->Px,chi->Py,chi->Pz,chi->E);
		Pachi.SetXYZT(achi->Px,achi->Py,achi->Pz,achi->E);
		
		//To get a factor 2 statistics, use both chi and achi
		int status1=chi->Status;
		int status2=achi->Status;
		
		
		hEnergy0->Fill(Pchi.E()); hEnergy0->Fill(Pachi.E());
		hAngle0->Fill(Pchi.Theta()*TMath::RadToDeg());hAngle0->Fill(Pachi.Theta()*TMath::RadToDeg());
		hSpace0->Fill(Pchi.Px()/Pchi.Pz(),Pchi.Py()/Pchi.Pz()); hSpace0->Fill(Pachi.Px()/Pachi.Pz(),Pachi.Py()/Pachi.Pz());
		hSpace1->Fill(10*Pchi.Px()/Pchi.Pz(),10*Pchi.Py()/Pchi.Pz()); hSpace1->Fill(10*Pachi.Px()/Pachi.Pz(),10*Pachi.Py()/Pachi.Pz());
		/*Here compute the acceptance*/
		//2 m
		x=(Pchi.Px()/Pchi.Pz())*200;
		y=(Pchi.Py()/Pchi.Pz())*200;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc2a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc2b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc2c++;		
		x=(Pachi.Px()/Pachi.Pz())*200;
		y=(Pachi.Py()/Pachi.Pz())*200;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc2a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc2b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc2c++;
		
		//5 m
		x=(Pchi.Px()/Pchi.Pz())*500;
		y=(Pchi.Py()/Pchi.Pz())*500;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc5a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc5b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc5c++;		
		x=(Pachi.Px()/Pachi.Pz())*500;
		y=(Pachi.Py()/Pachi.Pz())*500;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc5a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc5b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc5c++;
		
		//10 m
		x=(Pchi.Px()/Pchi.Pz())*1000;
		y=(Pchi.Py()/Pchi.Pz())*1000;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc10a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc10b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc10c++;
		x=(Pachi.Px()/Pachi.Pz())*1000;
		y=(Pachi.Py()/Pachi.Pz())*1000;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc10a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc10b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc10c++;
		
		//15 m
		x=(Pchi.Px()/Pchi.Pz())*1500;
		y=(Pchi.Py()/Pchi.Pz())*1500;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc15a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc15b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc15c++;	
		x=(Pachi.Px()/Pachi.Pz())*1500;
		y=(Pachi.Py()/Pachi.Pz())*1500;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc15a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc15b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc15c++;
  		
		//20 m
		x=(Pchi.Px()/Pchi.Pz())*2000;
		y=(Pchi.Py()/Pchi.Pz())*2000;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc20a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc20b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc20c++;	
		x=(Pachi.Px()/Pachi.Pz())*2000;
		y=(Pachi.Py()/Pachi.Pz())*2000;
		if ((fabs(x)<10)&&(fabs(y)<10))	 Nacc20a++;
		if ((fabs(x)<25)&&(fabs(y)<25))	 Nacc20b++;
		if ((fabs(x)<50)&&(fabs(y)<50))	 Nacc20c++;
		
		
		if ((status1>=1)&&(status2==0)) Pchi_in=&Pchi;
		else if ((status1==0)&&(status2>=1)) Pchi_in=&Pachi;
		else if ((status1==0)&&(status2==0)) continue;
		else if ((status1>=1)&&(status2>=1)){
			
			if (status1==11) Pchi_in=&Pchi;
			else if (status2==11) Pchi_in=&Pachi;
		}
		
		{
			hEnergy1->Fill(Pchi_in->E());
			hAngle1->Fill(Pchi_in->Theta()*TMath::RadToDeg());
			hVin0->Fill(vin.X(),vin.Y());
			hVout0->Fill(vout.X(),vout.Y());
			hL0->Fill(L/100.);
			
			
			hLvsE0->Fill(Pchi_in->E(),L/100.);
			hLvsTheta0->Fill(Pchi_in->Theta()*TMath::RadToDeg(),L/100.);
			
		}
		
		
		
	}
	
  //Draw
	TCanvas *c1=new TCanvas("c1","c1");
	c1->Divide(4,4);
	
	c1->cd(1);
	hEnergy0->Draw();
	hEnergy1->SetLineColor(2);
	hEnergy1->Draw("SAME");
	
	c1->cd(2);
	hAngle0->Draw();
	hAngle1->SetLineColor(2);
	hAngle1->Draw("SAME");
	
	c1->cd(3);
	hSpace0->Draw("colz");
	
	c1->cd(4);
	hSpace1->Draw("colz");
	
	c1->cd(5);
	hVin0->Draw("colz");
	
	c1->cd(6);
	hVout0->Draw("colz");
	
	c1->cd(7)->SetLogy();
	hL0->Draw();
	hL1->SetLineColor(2);
	hL2->SetLineColor(3);
	hL1->Draw("SAME");
	hL2->Draw("SAME");
	
	c1->cd(8);
	hLvsE0->Draw("colz");
	
	c1->cd(9);
	hLvsTheta0->Draw("colz");
	
	//Divide the acceptance counts
	Nacc2a/=(2*N);
	Nacc2b/=(2*N);
	Nacc2c/=(2*N);
	
	Nacc5a/=(2*N);
	Nacc5b/=(2*N);
	Nacc5c/=(2*N);
	
	Nacc10a/=(2*N);
	Nacc10b/=(2*N);
	Nacc10c/=(2*N);
	
	Nacc15a/=(2*N);
	Nacc15b/=(2*N);
	Nacc15c/=(2*N);
	
	Nacc20a/=(2*N);
	Nacc20b/=(2*N);
	Nacc20c/=(2*N);
	
	cout<<"Acceptance: "<<endl;
	cout<<Nacc2a<<" "<<Nacc2b<<" "<<Nacc2c<<endl;
	cout<<Nacc5a<<" "<<Nacc5b<<" "<<Nacc5c<<endl;
	cout<<Nacc10a<<" "<<Nacc10b<<" "<<Nacc10c<<endl;
	cout<<Nacc15a<<" "<<Nacc15b<<" "<<Nacc15c<<endl;
	cout<<Nacc20a<<" "<<Nacc20b<<" "<<Nacc20c<<endl;	
}


