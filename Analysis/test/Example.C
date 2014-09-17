//------------------------------------------------------------------------------

struct MyPlots
{
  TH1 *fTrackDeltaPT;
  TH1 *fJetPT[2];
  TH1 *fMissingET;
  TH1 *fElectronPT;
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, MyPlots *plots)
{
  THStack *stack;
  TLegend *legend;
  TPaveText *comment;

  // book 2 histograms for PT of 1st and 2nd leading jets

  plots->fJetPT[0] = result->AddHist1D("jet_pt_0", "leading jet P_{T}",
                                       "jet P_{T}, GeV/c", "number of jets",
                                       50, 0.0, 100.0);

  plots->fJetPT[1] = result->AddHist1D("jet_pt_1", "2nd leading jet P_{T}",
                                       "jet P_{T}, GeV/c", "number of jets",
                                       50, 0.0, 100.0);

  plots->fJetPT[0]->SetLineColor(kRed);
  plots->fJetPT[1]->SetLineColor(kBlue);

  // book 1 stack of 2 histograms

  stack = result->AddHistStack("jet_pt_all", "1st and 2nd jets P_{T}");
  stack->Add(plots->fJetPT[0]);
  stack->Add(plots->fJetPT[1]);

  // book legend for stack of 2 histograms

  legend = result->AddLegend(0.72, 0.86, 0.98, 0.98);
  legend->AddEntry(plots->fJetPT[0], "leading jet", "l");
  legend->AddEntry(plots->fJetPT[1], "second jet", "l");

  // attach legend to stack (legend will be printed over stack in .eps file)

  result->Attach(stack, legend);

  // book more histograms

  plots->fMissingET = result->AddHist1D("missing_et", "Missing E_{T}",
                                        "Missing E_{T}, GeV",
                                        "number of events",
                                        60, 0.0, 30.0);

  plots->fElectronPT = result->AddHist1D("electron_pt", "electron P_{T}",
                                         "electron P_{T}, GeV/c",
                                         "number of electrons",
                                         50, 0.0, 100.0);

  // book general comment

  comment = result->AddComment(0.54, 0.72, 0.98, 0.98);
  comment->AddText("demonstration plot");
  comment->AddText("produced by Example.C");

  // attach comment to single histograms

  result->Attach(plots->fJetPT[0], comment);
  result->Attach(plots->fJetPT[1], comment);
  result->Attach(plots->fMissingET, comment);
  result->Attach(plots->fElectronPT, comment);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
{
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  TRootJet *jet[2];
  TRootMissingET *met;
  TRootElectron *electron;

  TIter itElectron(branchElectron);

  Long64_t entry;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Analyse two leading jets
    if(branchJet->GetEntriesFast() >= 2)
    {
      jet[0] = (TRootJet*) branchJet->At(0);
      jet[1] = (TRootJet*) branchJet->At(1);

      plots->fJetPT[0]->Fill(jet[0]->PT);
      plots->fJetPT[1]->Fill(jet[1]->PT);
    }

    // Analyse missing ET
    if(branchMissingET->GetEntriesFast() > 0)
    {
      met = (TRootMissingET*) branchMissingET->At(0);
      plots->fMissingET->Fill(met->MET);
    }

    // Loop over all electrons in event
    itElectron.Reset();
    while((electron = (TRootElectron*) itElectron.Next()))
    {
      plots->fElectronPT->Fill(electron->PT);
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  result->Print();
}

//------------------------------------------------------------------------------

void Example(const char *inputFileList)
{
  TChain *chain = new TChain("LHCO");

  if(!FillChain(chain, inputFileList)) return;

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  MyPlots *plots = new MyPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------

