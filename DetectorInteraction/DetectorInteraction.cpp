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
//---------------------------------------------------------------------------

double AnalyseParticles(LHEF::Reader *reader) {
	LHEF::HEPEUP &hepeup = reader->hepeup;
	const LHEF::HEPRUP &heprup = reader->heprup;
	long PID;
	Int_t particle, n_inside;

	Double_t signPz, cosTheta,M;
	Double_t w,L,ndet,sigma;
	Double_t xthetaxfmin, xthetaxfmax, xthetayfmin, xthetayfmax;
	TLorentzVector chi, recoil_chi, recoil_elastic;
	TVector3 vin, vhit,vout, fiducialV;
	vector<double> tmp_v;
	vector<int> ii_inside;
	n_inside = 0;
	xthetaxfmin = xthetayfmin = 0;
	xthetaxfmax = (heprup.lx / (2 * heprup.ldet)); //ok without ATAN for this check
	xthetayfmax = (heprup.ly / (2 * heprup.ldet)); //ok without ATAN for this check

	//set the fiducial volume "box", with respect to the z axis.
	fiducialV.SetXYZ(heprup.lx, heprup.ly, heprup.lz);
	//init vhit
	vhit.SetXYZ(0.,0.,0.);
	//init vin
	vin.SetXYZ(0.,0.,0.);
	//init vout
	vout.SetXYZ(0.,0.,0.);

	w=hepeup.XWGTUP; //this is the event weight, in pbarn, as given by Madgraph. --> Cross section is the sum over the events of the event weight
	ndet=heprup.NDET;
	for (particle = 0; particle < hepeup.NUP; ++particle) {

		PID = hepeup.IDUP[particle];
		if ((PID != -611) && (PID != 611)){
			continue; //We are not interested in other particles. Go on
		}
		M = hepeup.PUP[particle][4];
		chi.SetPxPyPzE(hepeup.PUP[particle][0], hepeup.PUP[particle][1],
				hepeup.PUP[particle][2], hepeup.PUP[particle][3]);
		cosTheta = TMath::Abs(chi.CosTheta());
		signPz = (chi.Pz() >= 0.0) ? 1.0 : -1.0;
		/* I need now to apply the fiducial cuts.
		 I need to find which of the chis are inside (if more than one), and use that for the interaction.
		 If both are inside, I take only one.*/
		if ((fabs(chi.Px() / chi.Pz()) < xthetaxfmax)&& (fabs(chi.Py() / chi.Pz()) < xthetayfmax)){
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
		chi.SetPxPyPzE(hepeup.PUP[ii_inside.at(0)][0], hepeup.PUP[ii_inside.at(0)][1],
		               hepeup.PUP[ii_inside.at(0)][2], hepeup.PUP[ii_inside.at(0)][3]);
		
		vin.SetXYZ((chi.Px() / chi.Pz()) * heprup.ldet,
				(chi.Py() / chi.Pz()) * heprup.ldet, heprup.ldet); //the chi hit position in the fiducial volume front-face
		
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
			sigma=m_utils->doElasticRecoil(chi, recoil_elastic, recoil_chi,heprup.procid);		
			msigma+=sigma;
	
			if (sigma==0){
				w=0;
				hepeup.ISTUP[ii_inside.at(0)]=-11; //mark this negative
				break;
			}
			else{
				L=m_utils->findInteractionPoint(chi, fiducialV, vin,vout, vhit);
				L=L*100; //since the above returns it in m;
				mL+=L;
				w=w*L*heprup.NDET*sigma;
				//w=sigma;
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
	
	

	/*use the eventComments for the vertex location (in m, in the form x y z)
	use the eventComments also to report the "production weight * interaction probability * dump luminosity".
	In this way, the output file is self-consistent. To have the number of events per EOT, it is sufficient to sum these weights.
	*/
	w=w*n_inside*1E-36;//by multiplying per n_inside, automatically I correct for the fact I have two chis, potentially both in the detector.
	//the factor 1E-36 is the conversion pbarn ---> cm2
	w*=heprup.NDUMP*heprup.LDUMP;
	
	reader->eventComments  = Form("IN: %f %f %f \n",vin.X(),vin.Y(),vin.Z());
	reader->eventComments += Form("OUT: %f %f %f \n", vout.X(), vout.Y(), vout.Z());
	reader->eventComments += Form("HIT: %f %f %f \n", vhit.X(), vhit.Y(), vhit.Z());
	reader->eventComments += Form("L: %f \n",L);
	reader->eventComments += Form("W: %e",w);
	return w;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
	char appName[] = "BDX-DetectorInteraction";

	if (argc != 3) {
		cout << " Usage: " << appName << " input_file" << " output_file"
				<< endl;
		cout << " input_file - input file in LHEF format," << endl;
		cout << " output_file - output file in LHEF format." << endl;
		cout << " output_file is also used to write the number of events per EOT, as output_file.txt"<<endl;
		return 1;
	}

	gROOT->SetBatch();

	int appargc = 1;
	char *appargv[] = { appName };
	double W,L;
	double thisW;
	int Nin=0;
	// Open a stream connected to an event file:
	ifstream inputFileStream(argv[1]);
	ofstream outputFileStream(argv[2]);
	
	string outputFileEOTname(argv[2]);
	outputFileEOTname+=".summary.txt";
	ofstream outputFileEOT(outputFileEOTname.c_str());
	
	// Create the Reader object:
	LHEF::Reader *inputReader = new LHEF::Reader(inputFileStream); //this triggers also "init"
	// Create the Writer object:
	LHEF::Writer *outputWriter = new LHEF::Writer(outputFileStream);
	outputWriter->heprup = inputReader->heprup;
	outputWriter->headerStream.str(inputReader->headerBlock);
	outputWriter->init();

	//Instead of having, for each process, a single function with a free parameter (the incident chi energy),
	//it is much better to have 100 of them, each for a certain energy: otherwise, the integral gets calculated for each event!
	cout << " ** Preparing the cross-section functions and the KinUtils class ** " << endl;

	//note that alphaEM is saved in the event header (although we are not going to touch this at all!).
	//Therefore, I will set it in the first event in the  loop.
	m_utils=new KinUtils(inputReader->heprup.EBEAM,inputReader->heprup.FMASS,inputReader->heprup.APMASS,0,inputReader->heprup.EPSILON,inputReader->heprup.ALPHAD,inputReader->heprup.eTHR,inputReader->heprup.pTHR,inputReader->heprup.SEED);


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
			if (inputReader->heprup.procid)
				thisW=AnalyseParticles(inputReader); //this also returns the "corrected" event weight (production weight * interaction probability * dump luminosity)
				W+=thisW;
				if (thisW>0){
					Nin++;
				}
				outputWriter->hepeup = inputReader->hepeup;
				outputWriter->eventStream.str(inputReader->eventComments);
				outputWriter->writeEvent();
				progressBar.Update(entry);
				++entry;
		}
		progressBar.Finish();
	}
	cout<< " Events per EOT: "<<W<<endl;
	cout << "** Exiting..." << endl;
	
	cout<< " mean sigma: "<<msigma/Nin<<endl;
	cout<< " mean L: "<<mL/Nin<<endl;
	outputFileEOT<<W<<endl;
	outputFileEOT.close();
	
	delete inputReader;
	delete outputWriter;
	

}

