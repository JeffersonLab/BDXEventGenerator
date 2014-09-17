#include <iostream>

#include "TApplication.h"
#include "TLorentzVector.h"

#include "TFile.h"

#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

using namespace std;

// generated particle list

const int nmxhep = 4000;

typedef struct
{
  int nevhep;              // event number
  int nhep;                // number of entries in record
  int isthep[nmxhep];      // status code
  int idhep[nmxhep];       // particle ID (PDG standard)
  int jmohep[nmxhep][2];   // index to first and second particle mothers
  int jdahep[nmxhep][2];   // index to first and last daughter particles
  double phep[nmxhep][5];  // 4-vector and mass
  double vhep[nmxhep][4];  // (x,y,z) of production, and production time (mm/c)
} hepevtF77;

extern hepevtF77 hepevt_;

// PGS track list

const int ntrkmx = 500;

typedef struct
{
  int numtrk;              // number of tracks
  int dumtrk;              // number of tracks
  int indtrk[ntrkmx];      // index to HEPEVT particle
  double ptrk[ntrkmx][3];  // track 3-vector
  double qtrk[ntrkmx];     // track charge
} pgstrkF77;

extern pgstrkF77 pgstrk_;

// PGS calorimeter tower arrays

const int nphimax = 600;
const int netamax = 600;

typedef struct
{
  double ecal[nphimax][netamax];  // electromagnetic energy in each tower
  double hcal[nphimax][netamax];  // hadronic energy in each tower
  double met_cal;                 // calorimeter missing ET
  double phi_met_cal;             // calorimeter missing ET phi
} pgscalF77;

extern pgscalF77 pgscal_;


// PGS calorimeter cluster list

const int nclumx = 50;

typedef struct
{
  int cclu[nphimax][netamax];  // map of cluster indices
  int numclu;                  // number of clusters in list
  int dumclu;                  // number of clusters in list
  double pclu[nclumx][5];      // cluster 4 vector and mass
  int etaclu[nclumx];          // cluster seed tower eta
  int phiclu[nclumx];          // cluster seed tower phi
  double emclu[nclumx];        // cluster electromagnetic energy
  double ehclu[nclumx];        // cluster hadronic energy
  double efclu[nclumx];        // cluster electromagnetic fraction
  double widclu[nclumx];       // cluster width sqrt(deta**2+dphi**2)
  int mulclu[nclumx];          // cluster tower multiplicity
} pgscluF77;

extern pgscluF77 pgsclu_;

// PGS trigger object list

const int ntrgmx = 500;

typedef struct
{
  int numtrg;                 // number of trigger objects
  int dumtrg;                 // number of trigger objects
  int indtrg[ntrgmx];         // index to HEPEVT particle (where relevant)
  int typtrg[ntrgmx];         // reconstructed type:  0=photon
                              //                      1=electron
                              //                      2=muon
                              //                      3=tau (hadronic)
                              //                      4=jet
                              //                      5=detached vertex
                              //                      6=MET
  double vectrg[ntrgmx][10];  // trigger object vector: 1 = eta
                              //                        2 = phi
                              //                        3 = ET of cluster
                              //                        4 = cluster #
                              //                        5 = pt of track (if any)
                              //                        6 = track #
} pgstrgF77;

extern pgstrgF77 pgstrg_;

// PGS reconstructed object list

const int nobjmx = 500;

typedef struct
{
  int numobj;                 // number of reconstructed objects
  int dumobj;                 // number of reconstructed objects
  int indobj[nobjmx];         // index to HEPEVT particle (where relevant)
  int typobj[nobjmx];         // reconstructed type:  0 = photon
                              //                      1 = electron
                              //                      2 = muon
                              //                      3 = tau (hadronic)
                              //                      4 = jet
                              //                      5 = heavy charged
  double pobj[nobjmx][4];     // four vector of reconstructed object
  double qobj[nobjmx];        // charge of reconstructed object
  double vecobj[nobjmx][10];  // interesting object quantities
} pgsrecF77;

extern pgsrecF77 pgsrec_;

//           --------------------------
//           table of vecobj quantities
//           --------------------------
//
//   -------------------------------------------------------------------------------------
//    type            1           2          3        4        5        6        7
//   -------------------------------------------------------------------------------------
//   0  photon     EM energy  HAD energy  track E   N(trk)   width      -        -
//   1  electron    "   "      "     "       "        "        -        -        -
//   2  muon        "   "      "     "       "        "     trk iso E   -        -
//   3  tau         "   "      "     "       "        "      width    mtau     ptmax
//   4  jet         "   "      "     "       "        "        "      flavor   c,b tags ->
//   -------------------------------------------------------------------------------------
//
// b, c tagging: vecobj(7,iobj) non-zero if charm tag (jet prob. alg.)
//               vecobj(8,iobj) non-zero if b     tag ( "    "    "  )
//               vecobj(9,iobj) non-zero if b     tag (impact method)
//
//   --> all algorithms include rates for gluon, uds, c and b jets
//

static TFile *outputFile;
static ExRootTreeWriter *treeWriter;

static ExRootTreeBranch *branchGenParticle;
static ExRootTreeBranch *branchTrack;
static ExRootTreeBranch *branchCalTower;
static ExRootTreeBranch *branchMissingET;
static ExRootTreeBranch *branchCalCluster;
static ExRootTreeBranch *branchPhoton;
static ExRootTreeBranch *branchElectron;
static ExRootTreeBranch *branchMuon;
static ExRootTreeBranch *branchTau;
static ExRootTreeBranch *branchJet;
static ExRootTreeBranch *branchHeavy;

static void analyse_particle(Int_t number, ExRootTreeBranch *branch);
static void analyse_track(Int_t number, ExRootTreeBranch *branch);
static void analyse_tower(Int_t eta, Int_t phi, ExRootTreeBranch *branch);
static void analyse_met(ExRootTreeBranch *branch);
static void analyse_cluster(Int_t number, ExRootTreeBranch *branch);
static void analyse_photon(Int_t number, ExRootTreeBranch *branch);
static void analyse_electron(Int_t number, ExRootTreeBranch *branch);
static void analyse_muon(Int_t number, ExRootTreeBranch *branch);
static void analyse_tau(Int_t number, ExRootTreeBranch *branch);
static void analyse_jet(Int_t number, ExRootTreeBranch *branch);
static void analyse_heavy(Int_t number, ExRootTreeBranch *branch);

extern "C"
{
  void test_cpp__(int *int_var)
  {
    cout << "int_var = " << *int_var << endl;

    cout << "nevhep = " << hepevt_.nevhep << endl;
    cout << "jmohep(1,2) = " << hepevt_.jmohep[1][0] << endl;
    cout << "phep(1,1) = " << hepevt_.phep[0][0] << endl;

    cout << "numtrk = " << pgstrk_.numtrk << endl;
    cout << "ptrk(2,1) = " << pgstrk_.ptrk[0][1] << endl;

    cout << "hcal(2,1) = " << pgscal_.hcal[0][1] << endl;
    cout << "met_cal = " << pgscal_.met_cal << endl;

    cout << "pclu(2,1) = " << pgsclu_.pclu[0][1] << endl;
    cout << "mulclu(2) = " << pgsclu_.mulclu[1] << endl;

    cout << "vectrg(9,2) = " << pgstrg_.vectrg[1][8] << endl;
    cout << "indtrg(2) = " << pgstrg_.indtrg[1] << endl;

    cout << "vecobj(10,3) = " << pgsrec_.vecobj[2][9] << endl;
    cout << "indobj(3) = " << pgsrec_.indobj[2] << endl;
  }

//------------------------------------------------------------------------------

  void pgs2root_ini__()
  {
    int appargc = 2;
    char *appargv[] = {"pgs2root", "-b"};
    TApplication app("pgs2root", &appargc, appargv);
  
    TString outputFileName("pgs.root");
    TString treeName("PGS");
  
    outputFile = TFile::Open(outputFileName, "RECREATE");
    treeWriter = new ExRootTreeWriter(outputFile, treeName);

    // generated particles from HEPEVT
    branchGenParticle = treeWriter->NewBranch("GenParticle", TRootGenParticle::Class());
    // reconstructed tracks
    branchTrack = treeWriter->NewBranch("Track", TRootTrack::Class());
    // reconstructed calorimeter towers
    branchCalTower = treeWriter->NewBranch("CalTower", TRootCalTower::Class());
    // missing transverse energy
    branchMissingET = treeWriter->NewBranch("MissingET", TRootMissingET::Class());
    // reconstructed calorimeter clusters for jets and tau leptons
    branchCalCluster = treeWriter->NewBranch("CalCluster", TRootCalCluster::Class());
    // reconstructed photons
    branchPhoton = treeWriter->NewBranch("Photon", TRootPhoton::Class());
    // reconstructed electrons
    branchElectron = treeWriter->NewBranch("Electron", TRootElectron::Class());
    // reconstructed muons
    branchMuon = treeWriter->NewBranch("Muon", TRootMuon::Class());
    // reconstructed tau leptons
    branchTau = treeWriter->NewBranch("Tau", TRootTau::Class());
    // reconstructed jets
    branchJet = treeWriter->NewBranch("Jet", TRootJet::Class());
    // reconstructed heavy particles
    branchHeavy = treeWriter->NewBranch("Heavy", TRootHeavy::Class());
  }

//---------------------------------------------------------------------------

  void pgs2root_evt__()
  {
    Int_t particle, track, eta, phi, cluster, object;
  
    treeWriter->Clear();
  
    for(particle = 0; particle < hepevt_.nhep; ++particle)
    {
      analyse_particle(particle, branchGenParticle);
    }
  
    for(track = 0; track < pgstrk_.numtrk; ++track)
    {
      analyse_track(track, branchTrack);
    }
  
    for(eta = 0; eta < netamax; ++eta)
    {
      for(phi = 0; phi < nphimax; ++phi)
      {
        if(pgscal_.ecal[phi][eta] == 0.0 &&
           pgscal_.hcal[phi][eta] == 0.0) continue;

        analyse_tower(eta, phi, branchCalTower);
      }
    }

    analyse_met(branchMissingET);
    
    for(cluster = 0; cluster < pgsclu_.numclu; ++ cluster)
    {
      analyse_cluster(cluster, branchCalCluster);
    }
  
    for(object = 0; object < pgsrec_.numobj; ++object)
    {
      switch(pgsrec_.typobj[object])
      {
        case 0: analyse_photon(object, branchPhoton); break;
        case 1: analyse_electron(object, branchElectron); break;
        case 2: analyse_muon(object, branchMuon); break;
        case 3: analyse_tau(object, branchTau); break;
        case 4: analyse_jet(object, branchJet); break;
        case 5: analyse_heavy(object, branchHeavy); break;
      }
    }
    
    treeWriter->Fill();
  }

//---------------------------------------------------------------------------

  void pgs2root_end__()
  {
    treeWriter->Write();
    
    delete treeWriter;
    delete outputFile;
  }

} // extern "C"

//---------------------------------------------------------------------------

static void analyse_particle(Int_t number, ExRootTreeBranch *branch)
{
  TRootGenParticle *entry;

  Double_t signEta;

  entry = static_cast<TRootGenParticle*>(branch->NewEntry());

  entry->PID = hepevt_.idhep[number];
  entry->Status = hepevt_.isthep[number];
  entry->M1 = hepevt_.jmohep[number][0] - 1;
  entry->M2 = hepevt_.jmohep[number][1] - 1;
  entry->D1 = hepevt_.jdahep[number][0] - 1;
  entry->D2 = hepevt_.jdahep[number][1] - 1;

  entry->E = hepevt_.phep[number][3];
  entry->Px = hepevt_.phep[number][0];
  entry->Py = hepevt_.phep[number][1];
  entry->Pz = hepevt_.phep[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->T = hepevt_.vhep[number][3];
  entry->X = hepevt_.vhep[number][0];
  entry->Y = hepevt_.vhep[number][1];
  entry->Z = hepevt_.vhep[number][2];
}

//---------------------------------------------------------------------------

static void analyse_track(Int_t number, ExRootTreeBranch *branch)
{
  TRootTrack *entry;

  Double_t signEta;

  entry = static_cast<TRootTrack*>(branch->NewEntry());

  entry->Px = pgstrk_.ptrk[number][0];
  entry->Py = pgstrk_.ptrk[number][1];
  entry->Pz = pgstrk_.ptrk[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->Charge = pgstrk_.qtrk[number];

  entry->ParticleIndex = pgstrk_.indtrk[number] - 1;  
}

//---------------------------------------------------------------------------

static void analyse_tower(Int_t eta, Int_t phi, ExRootTreeBranch *branch)
{
  TRootCalTower *entry;

  entry = static_cast<TRootCalTower*>(branch->NewEntry());

  entry->Eem = pgscal_.ecal[phi][eta];
  entry->Ehad = pgscal_.hcal[phi][eta];
  entry->E = entry->Eem + entry->Ehad;
}

//---------------------------------------------------------------------------

static void analyse_met(ExRootTreeBranch *branch)
{
  TRootMissingET *entry;

  entry = static_cast<TRootMissingET*>(branch->NewEntry());

  entry->MET = pgscal_.met_cal;
  entry->Phi = pgscal_.phi_met_cal;
}

//---------------------------------------------------------------------------

static void analyse_cluster(Int_t number, ExRootTreeBranch *branch)
{
  TRootCalCluster *entry;

  entry = static_cast<TRootCalCluster*>(branch->NewEntry());

  entry->E = pgsclu_.pclu[number][3];
  entry->Px = pgsclu_.pclu[number][0];
  entry->Py = pgsclu_.pclu[number][1];
  entry->Pz = pgsclu_.pclu[number][2];

  entry->Eta = pgsclu_.etaclu[number];
  entry->Phi = pgsclu_.phiclu[number];

  entry->Eem = pgsclu_.emclu[number];
  entry->Ehad = pgsclu_.ehclu[number];
  entry->EemOverEtot = pgsclu_.efclu[number];

  entry->Ntwr = pgsclu_.mulclu[number];
}

//---------------------------------------------------------------------------

static void analyse_photon(Int_t number, ExRootTreeBranch *branch)
{
  TRootPhoton *entry;

  Double_t signEta;
  
  entry = static_cast<TRootPhoton*>(branch->NewEntry());

  entry->E = pgsrec_.pobj[number][3];
  entry->Px = pgsrec_.pobj[number][0];
  entry->Py = pgsrec_.pobj[number][1];
  entry->Pz = pgsrec_.pobj[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->Eem = pgsrec_.vecobj[number][0];
  entry->Ehad = pgsrec_.vecobj[number][1];
  entry->PTtrk = pgsrec_.vecobj[number][2];

  entry->Niso = pgsrec_.vecobj[number][3];

  entry->ET = pgsrec_.vecobj[number][5];
  entry->ETiso = pgsrec_.vecobj[number][6];
  entry->PTiso = pgsrec_.vecobj[number][7];

  entry->EhadOverEem = pgsrec_.vecobj[number][8];
  entry->EemOverPtrk = pgsrec_.vecobj[number][9];

  entry->ParticleIndex = pgsrec_.indobj[number] - 1;  
}

//---------------------------------------------------------------------------

static void analyse_electron(Int_t number, ExRootTreeBranch *branch)
{
  TRootElectron *entry;

  Double_t signEta;

  entry = static_cast<TRootElectron*>(branch->NewEntry());

  entry->E = pgsrec_.pobj[number][3];
  entry->Px = pgsrec_.pobj[number][0];
  entry->Py = pgsrec_.pobj[number][1];
  entry->Pz = pgsrec_.pobj[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->Charge = pgsrec_.qobj[number];

  entry->Eem = pgsrec_.vecobj[number][0];
  entry->Ehad = pgsrec_.vecobj[number][1];
  entry->PTtrk = pgsrec_.vecobj[number][2];

  entry->Niso = pgsrec_.vecobj[number][3];

  entry->ET = pgsrec_.vecobj[number][5];
  entry->ETiso = pgsrec_.vecobj[number][6];
  entry->PTisoMinusPTtrk = pgsrec_.vecobj[number][7];

  entry->EhadOverEem = pgsrec_.vecobj[number][8];
  entry->EemOverPtrk = pgsrec_.vecobj[number][9];

  entry->ParticleIndex = pgsrec_.indobj[number] - 1;
}

//---------------------------------------------------------------------------

static void analyse_muon(Int_t number, ExRootTreeBranch *branch)
{
  TRootMuon *entry;

  Double_t signEta;

  entry = static_cast<TRootMuon*>(branch->NewEntry());

  entry->E = pgsrec_.pobj[number][3];
  entry->Px = pgsrec_.pobj[number][0];
  entry->Py = pgsrec_.pobj[number][1];
  entry->Pz = pgsrec_.pobj[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->Charge = pgsrec_.qobj[number];

  entry->Eem = pgsrec_.vecobj[number][0];
  entry->Ehad = pgsrec_.vecobj[number][1];
  entry->Ptrk = pgsrec_.vecobj[number][2];

  entry->Niso = pgsrec_.vecobj[number][3];

  entry->Etrk = pgsrec_.vecobj[number][4];
  entry->ETiso = pgsrec_.vecobj[number][5];
  entry->PTiso = pgsrec_.vecobj[number][6];

  entry->ParticleIndex = pgsrec_.indobj[number] - 1;
}

//---------------------------------------------------------------------------

static void analyse_tau(Int_t number, ExRootTreeBranch *branch)
{
  TRootTau *entry;

  Double_t signEta;

  entry = static_cast<TRootTau*>(branch->NewEntry());

  entry->E = pgsrec_.pobj[number][3];
  entry->Px = pgsrec_.pobj[number][0];
  entry->Py = pgsrec_.pobj[number][1];
  entry->Pz = pgsrec_.pobj[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->Charge = pgsrec_.qobj[number];

  entry->Eem = pgsrec_.vecobj[number][0];
  entry->Ehad = pgsrec_.vecobj[number][1];
  entry->Etrk = pgsrec_.vecobj[number][2];

  entry->Ntrk = pgsrec_.vecobj[number][3];

  entry->Width = pgsrec_.vecobj[number][4];
  entry->Ecut = pgsrec_.vecobj[number][5];
  entry->PTmax = pgsrec_.vecobj[number][6];

  entry->SeedTrackIndex = pgsrec_.vecobj[number][7] - 1;
  entry->ClusterIndex = pgsrec_.indobj[number] - 1;
}

//---------------------------------------------------------------------------

static void analyse_jet(Int_t number, ExRootTreeBranch *branch)
{
  TRootJet *entry;

  Double_t signEta;

  entry = static_cast<TRootJet*>(branch->NewEntry());

  entry->E = pgsrec_.pobj[number][3];
  entry->Px = pgsrec_.pobj[number][0];
  entry->Py = pgsrec_.pobj[number][1];
  entry->Pz = pgsrec_.pobj[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();

  entry->Charge = pgsrec_.qobj[number];

  entry->Eem = pgsrec_.vecobj[number][0];
  entry->Ehad = pgsrec_.vecobj[number][1];
  entry->Etrk = pgsrec_.vecobj[number][2];

  entry->Ntrk = pgsrec_.vecobj[number][3];

  entry->Width = pgsrec_.vecobj[number][4];

  entry->Type = pgsrec_.vecobj[number][5]; // type: 21=g, 1=d, 2=u, 3=s, 4=c, 5=b

  entry->CTag = pgsrec_.vecobj[number][6];
  entry->BTagVtx = pgsrec_.vecobj[number][7];
  entry->BTagImp = pgsrec_.vecobj[number][8];

  entry->ClusterIndex = pgsrec_.indobj[number] - 1;  
}

//---------------------------------------------------------------------------

static void analyse_heavy(Int_t number, ExRootTreeBranch *branch)
{
  TRootHeavy *entry;

  Double_t signEta;

  entry = static_cast<TRootHeavy*>(branch->NewEntry());

  entry->E = pgsrec_.pobj[number][3];
  entry->Px = pgsrec_.pobj[number][0];
  entry->Py = pgsrec_.pobj[number][1];
  entry->Pz = pgsrec_.pobj[number][2];

  TVector3 vector(entry->Px, entry->Py, entry->Pz);

  entry->PT = vector.Perp();
  signEta = (entry->Pz >= 0.0) ? 1.0 : -1.0;
  entry->Eta = vector.CosTheta()*vector.CosTheta() == 1.0 ? signEta*999.9 : vector.Eta();
  entry->Phi = vector.Phi();
 
  entry->ParticleIndex = pgsrec_.indobj[number] - 1;
}

//---------------------------------------------------------------------------

