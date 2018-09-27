#include <iostream>
#include <map>
#include <cstdlib>      // std::rand, std::srand
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
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

#include "ExRootClasses.h"
#include "ExRootTreeWriter.h"
#include "ExRootTreeBranch.h"

#include "ExRootUtilities.h"
#include "ExRootProgressBar.h"

#include "KinUtils.h"


using namespace std;

KinUtils *m_utils;

double msigma=0;
double mL=0;

/*These are here for writeLund*/
TVector3 vin, vhit,vout, fiducialV,vMC;

//---------------------------------------------------------------------------

void writeLund(ofstream &ofile,LHEF::HEPEUP &data){
	int nparticles=data.IDUP.size();
	int pdg;
	int isGood=0;
	int status;
	int idParticle,idInitialChi,idFinalChi;
	double Px,Py,Pz,E,M,vx,vy,vz;

	int charge;
	/*We write to the lund file only the final state proton (or electron), the input chi, and the output chi.
	 * The latter twos have status==0
	 */
	/*1: preliminary check to see we have all the particles*/
	for (int ip=0;ip<nparticles;ip++){
		pdg=data.IDUP[ip];
		if ((pdg==92212)||(pdg==911)){
			isGood+=1;
			idParticle=ip;
			if (pdg==911) charge=-1;
			else if (pdg==92212) charge=1;

			vx=vhit.X()*100;
			vy=vhit.Y()*100;
			vz=vhit.Z()*100; /*in cm*/

		}
		else if (pdg==9611){
			isGood+=2;
			idFinalChi=ip;
		}
		else if ((pdg==611)|(pdg==-611)){
			status=data.ISTUP[ip];
			if (status==11){
				isGood+=4;
				idInitialChi=ip;
			}
		}
	}
	if (isGood!=7){
		return; /*For same reasons, in this events not all the particles where there??*/
	}
	else{
		ofile<<"3 0 0 0 0 0 0 0 0 0"<<endl; //lund header. 10 numbers, first is number of particles (the only relevant!)

		/*The initial state chi*/
		Px=data.PUP[idInitialChi][0];
		Py=data.PUP[idInitialChi][1];
		Pz=data.PUP[idInitialChi][2];
		E=data.PUP[idInitialChi][3];
		M=data.PUP[idInitialChi][4];
		ofile<<"1 0 0 611 0 0 "<<Px<<" "<<Py<<" "<<Pz<<" "<<E<<" "<<M<<" "<<vx<<" "<<vy<<" "<<vz<<endl;


		/*The final state chi*/
		Px=data.PUP[idFinalChi][0];
		Py=data.PUP[idFinalChi][1];
		Pz=data.PUP[idFinalChi][2];
		E=data.PUP[idFinalChi][3];
		M=data.PUP[idFinalChi][4];
		ofile<<"2 0 0 9611 0 0 "<<Px<<" "<<Py<<" "<<Pz<<" "<<E<<" "<<M<<" "<<vx<<" "<<vy<<" "<<vz<<endl;

		/*The final state SM particle - e- or P or nucleus*/
		Px=data.PUP[idParticle][0];
		Py=data.PUP[idParticle][1];
		Pz=data.PUP[idParticle][2];
		E=data.PUP[idParticle][3];
		M=data.PUP[idParticle][4];
		pdg=data.IDUP[idParticle];
		if (pdg==911) pdg=11;
		else if (pdg==92212) pdg=2212;
		else if (pdg==999) pdg=999; //the nucleus TODO. For now 999
		ofile<<"3 "<<charge<<" 1 "<<pdg<<" 0 0 "<<Px<<" "<<Py<<" "<<Pz<<" "<<E<<" "<<M<<" "<<vx<<" "<<vy<<" "<<vz<<endl;
	}
}

//Return: the production weight (first) and the interaction weight(second)
std::pair<double,double> AnalyseParticles(LHEF::Reader *reader) {


	LHEF::HEPEUP &hepeup = reader->hepeup;              /*This is a reference!*/
	const LHEF::HEPRUP &heprup = reader->heprup; 		/*This is a reference!*/
	long PID;
	Int_t particle, n_inside;

	Double_t signPz, cosTheta,sinTheta,cosPhi,sinPhi,M;
	Double_t w,L,ndet,sigma;
	Double_t xmin,xmax,ymin,ymax;
	Double_t x,y;
	TLorentzVector chi, recoil_chi, recoil_elastic;

	vector<double> tmp_v;
	vector<int> ii_inside;
	n_inside = 0;


	xmin=-heprup.lx/2+heprup.displacement;
	xmax=heprup.lx/2+heprup.displacement;
	ymin=-heprup.ly/2;
	ymax=heprup.ly/2;




	//set the fiducial volume "box", with respect to the z axis.
	fiducialV.SetXYZ(heprup.lx, heprup.ly, heprup.lz);
	//init vhit
	vhit.SetXYZ(0.,0.,0.);
	//init vin
	vin.SetXYZ(0.,0.,0.);
	//init vout
	vout.SetXYZ(0.,0.,0.);

	//Set the coordinate of the fiducial volume front face center in MC-Geant4 coordinates
	vMC.SetXYZ(heprup.MCcenterX,heprup.MCcenterY,heprup.MCcenterZ); //by definition this point, in the "MadGraph" reference, is at (0,0,ldet)


	w=hepeup.XWGTUP; //this is the event weight, in pbarn, as given by Madgraph. --> Cross section is the sum over the events of the event weight
	ndet=heprup.NDET;
	for (particle = 0; particle < hepeup.NUP; ++particle) {

		PID = hepeup.IDUP[particle];
		if ((PID != -611) && (PID != 611)){
			continue; //We are not interested in other particles. Go on
		}
		M = hepeup.PUP[particle][4];
		chi.SetPxPyPzE(hepeup.PUP[particle][0], hepeup.PUP[particle][1],hepeup.PUP[particle][2], hepeup.PUP[particle][3]);
		cosTheta = TMath::Abs(chi.CosTheta());
		sinTheta = sqrt(1-cosTheta*cosTheta);

		cosPhi = cos(chi.Phi());
		sinPhi = sin(chi.Phi());

		x=heprup.ldet*(sinTheta/cosTheta)*cosPhi;
		y=heprup.ldet*(sinTheta/cosTheta)*sinPhi;

	



		signPz = (chi.Pz() >= 0.0) ? 1.0 : -1.0;
		/* I need now to apply the fiducial cuts.
		 I need to find which of the chis are inside (if more than one), and use that for the interaction.
		 If both are inside, I take only one.*/

		if ((x>xmin)&&(x<xmax)&&(y>ymin)&&(y<ymax)){
			n_inside++;
			ii_inside.push_back(particle);
			hepeup.ISTUP[particle] = 1;

		}

		else { //for now, just mark this chi out of the fiducial volume with a status "0"
			hepeup.ISTUP[particle] = 0;
			continue;
		}
	}



	//now we know how many chis are inside (n_inside). Take 1 random.
	if (n_inside!=0){
		std::random_shuffle (ii_inside.begin(), ii_inside.end() );
		chi.SetPxPyPzE(hepeup.PUP[ii_inside.at(0)][0], hepeup.PUP[ii_inside.at(0)][1],hepeup.PUP[ii_inside.at(0)][2],hepeup.PUP[ii_inside.at(0)][3]);

		cosTheta = TMath::Abs(chi.CosTheta());
		sinTheta = sqrt(1-cosTheta*cosTheta);
		cosPhi = cos(chi.Phi());
		sinPhi = sin(chi.Phi());
		x=heprup.ldet*(sinTheta/cosTheta)*cosPhi;
		y=heprup.ldet*(sinTheta/cosTheta)*sinPhi;

		vin.SetXYZ(x,y,heprup.ldet); //the chi hit position in the fiducial volume front-face

		//mark this chi with status "11":
		hepeup.ISTUP[ii_inside.at(0)]=11;

		//From here, we have a chi within the fiducial volume.
		//See which interaction to consider
		switch (heprup.procid) {
		case Proc_nothing: //nothing to do
			w=0;
			break;
		case Proc_Pelastic: //proton elastic
		case Proc_Eelastic: //electron elastic
		case Proc_Nuclelastic: //nuclear elastic (coherent)
			sigma=m_utils->doElasticRecoil(chi, recoil_elastic, recoil_chi,heprup.procid);		
			msigma+=sigma;

			if (sigma==0){
				w=0;
				hepeup.ISTUP[ii_inside.at(0)]=-11; //mark this negative
				break;
			}
			else{

				/*No need to re-write displacement routine!*/
				vin.SetX(vin.X()-heprup.displacement);
				L=m_utils->findInteractionPoint(chi, fiducialV,vin,vout, vhit);

				vin.SetX(vin.X()+heprup.displacement);
				vout.SetX(vout.X()+heprup.displacement);
				vhit.SetX(vhit.X()+heprup.displacement);
				L=L*100; //since the above returns it in m;
				mL+=L;
				w=w*L*heprup.NDET*sigma;      /*Multiply the total cross-section * interaction length of this event * weight (pbarn) of this event*/

				//correct the particles positions in Geant4-MC format//
				vin.SetXYZ(vin.X()+vMC.X(),vin.Y()+vMC.Y(),vin.Z()+vMC.Z()-heprup.ldet);
				vout.SetXYZ(vout.X()+vMC.X(),vout.Y()+vMC.Y(),vout.Z()+vMC.Z()-heprup.ldet);
				vhit.SetXYZ(vhit.X()+vMC.X(),vhit.Y()+vMC.Y(),vhit.Z()+vMC.Z()-heprup.ldet);

				//add particles to hepeup
				//final state chi
				hepeup.IDUP.push_back(9611);
				hepeup.ISTUP.push_back(1);
				hepeup.MOTHUP.push_back(std::make_pair(particle + 1, particle + 1));
				hepeup.ICOLUP.push_back(std::make_pair(0, 0));
				tmp_v.clear();
				tmp_v.push_back(recoil_chi.Px());
				tmp_v.push_back(recoil_chi.Py());
				tmp_v.push_back(recoil_chi.Pz());
				tmp_v.push_back(recoil_chi.E());
				tmp_v.push_back(recoil_chi.M());
				hepeup.PUP.push_back(tmp_v);
				hepeup.VTIMUP.push_back(0);
				hepeup.SPINUP.push_back(0);
				//final state recoil
				(heprup.procid == Proc_Pelastic ? hepeup.IDUP.push_back(92212) : hepeup.IDUP.push_back(911));
				hepeup.ISTUP.push_back(1);
				hepeup.MOTHUP.push_back(std::make_pair(particle + 1, particle + 1));
				hepeup.ICOLUP.push_back(std::make_pair(0, 0));
				tmp_v.clear();
				tmp_v.push_back(recoil_elastic.Px());
				tmp_v.push_back(recoil_elastic.Py());
				tmp_v.push_back(recoil_elastic.Pz());
				tmp_v.push_back(recoil_elastic.E());
				tmp_v.push_back(recoil_elastic.M());
				hepeup.PUP.push_back(tmp_v);
				hepeup.VTIMUP.push_back(0);
				hepeup.SPINUP.push_back(0);
				hepeup.NUP += 2;
				break;
			}
		case Proc_Pinelastic: //inelastic
		case Proc_Einelastic:
			w=0;
			break;
		default:
			w=0;
			cout << "Error, interaction not recognized" << endl;
			break;
		}
	}
	else{ //ninside=0;
		w=0;
	}
	/*use the eventComments for the vertex location (in m, in the form x y z)
	use the eventComments also to report the "production weight * interaction probability * dump luminosity".
	In this way, the output file is self-consistent. To have the number of events per EOT, it is sufficient to sum these weights.
	 */
	w=w*n_inside*1E-36;//by multiplying per n_inside, automatically I correct for the fact I have two chis, potentially both in the detector.
	//the factor 1E-36 is the conversion pbarn ---> cm2
	w*=heprup.NDUMP*heprup.LDUMP;

	//also write the PRODUCTION cross section, multiplied by NDUMP and LDUMP, to have the number of total chi-chibar pairs produced by summing over this number
	sigma=hepeup.XWGTUP*heprup.NDUMP*heprup.LDUMP*1E-36; 	//the factor 1E-36 is the conversion pbarn ---> cm2

	reader->eventComments  = Form("IN: %f %f %f \n",vin.X(),vin.Y(),vin.Z());
	reader->eventComments += Form("OUT: %f %f %f \n", vout.X(), vout.Y(), vout.Z());
	reader->eventComments += Form("HIT: %f %f %f \n", vhit.X(), vhit.Y(), vhit.Z());
	reader->eventComments += Form("L: %f \n",L);
	reader->eventComments += Form("W: %e \n",w);
	reader->eventComments += Form("Sigma: %e",sigma);

	std::pair<double,double> ret=std::make_pair(sigma,w);
	return ret;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
	char appName[] = "DetectorInteraction.exe";

	if (argc != 3) {
		cout << " Usage: " << appName << " input_file" << " output_file_name"
				<< endl;
		cout << " input_file - input file (LHEF format)" << endl;
		cout << " output_file - output file name. Will create ouput_file.lhe and ouput_file.lund" << endl;
		cout << " output_file is also used to write the number of events per EOT, as output_file.txt"<<endl;
		return 1;
	}

	gROOT->SetBatch();

	int appargc = 1;
	char *appargv[] = { appName };
	double W,L,Sigma;
	double thisW,thisSigma;
	std::pair<double,double> retVal;
	int Nin=0;
	// Open a stream connected to an event file:
	string outputFileEOTName(argv[2]);
	outputFileEOTName+=".summary.txt";
	string outputFileLundName(argv[2]);
	outputFileLundName+=".lund";
	string outputFileLHEName(argv[2]);
	outputFileLHEName+=".lhe";

	cout<<"Out files: "<<endl;
	cout<<"LUND: "<<outputFileLundName<<endl;
	cout<<"LHE : "<<outputFileLHEName<<endl;
	cout<<"SUM : "<<outputFileEOTName<<endl;
	ifstream inputFile(argv[1]);

	ofstream outputFileLHE(outputFileLHEName.c_str());
	ofstream outputFileLund(outputFileLundName.c_str());
	ofstream outputFileEOT(outputFileEOTName.c_str());

	// Create the Reader object:
	LHEF::Reader *inputReader = new LHEF::Reader(inputFile); //this triggers also "init"
	// Create the Writer object:
	LHEF::Writer *outputWriter = new LHEF::Writer(outputFileLHE);
	outputWriter->heprup = inputReader->heprup;
	outputWriter->headerStream.str(inputReader->headerBlock);
	outputWriter->init();

	//Instead of having, for each process, a single function with a free parameter (the incident chi energy),
	//it is much better to have 100 of them, each for a certain energy: otherwise, the integral gets calculated for each event!
	cout << " ** Preparing the cross-section functions and the KinUtils class ** " << endl;

	//note that alphaEM is saved in the event header (although we are not going to touch this at all!).
	//Therefore, I will set it in the first event in the  loop.
	m_utils=new KinUtils(inputReader->heprup.EBEAM,inputReader->heprup.FMASS,inputReader->heprup.APMASS,0,inputReader->heprup.EPSILON,inputReader->heprup.ALPHAD,inputReader->heprup.eTHR,inputReader->heprup.pTHR,inputReader->heprup.pBINDING,inputReader->heprup.SEED);


	cout << "** Calculating number of events to process. Please wait..."
			<< endl;
	Long64_t allEntries = inputReader->getNumberOfEvents();
	cout << "** Input file contains " << allEntries << " events" << endl;
	W=0;
	if (allEntries > 0) {
		ExRootProgressBar progressBar(allEntries);

		// Loop over all events
		Long64_t entry = 0;
		while (inputReader->readEvent()) {

			if (entry==0){
				cout<<"First entry"<<endl;
				m_utils->setAlpha(inputReader->hepeup.AQEDUP); //set alpha_EM.
				gRandom->SetSeed(inputReader->heprup.SEED);
				std::srand(inputReader->heprup.SEED);

			}
			//This is the function that triggers the interaction in the fiducial volume.
			thisW=0;
			thisSigma=0;
			if (inputReader->heprup.procid){

				retVal=AnalyseParticles(inputReader); //this also returns the "corrected" event weight (production weight * interaction probability * dump luminosity)
				thisW=retVal.second;
				thisSigma=retVal.first;
			}
			W+=thisW;
			Sigma+=thisSigma;
			if (thisW>0){
				Nin++;
			}
			outputWriter->hepeup = inputReader->hepeup;
			outputWriter->eventStream.str(inputReader->eventComments);
			outputWriter->writeEvent();
			/*Now write the LUND block*/
			if (thisW>0) writeLund(outputFileLund,inputReader->hepeup);

			progressBar.Update(entry);
			++entry;

		}
		progressBar.Finish();
	}
	cout << "** Exiting..." << endl;
	cout<<" Original number of events in the file: "<<allEntries<<endl;
	cout<<" Number  of entries actually written in output: "<<Nin<<endl;
	cout<< " Events per EOT: "<<W<<endl;

	cout<< " mean sigma: "<<msigma/Nin<<endl;
	cout<< " mean L: "<<mL/Nin<<endl;

	outputFileEOT<<"Number of produced A' (chi pairs) per EOT: "<<endl;
	outputFileEOT<<Sigma<<endl;
	outputFileEOT<<"Number of interacting events per EOT: "<<endl;
	outputFileEOT<<W<<endl;
	outputFileEOT.close();

	delete inputReader;
	delete outputWriter;


}

