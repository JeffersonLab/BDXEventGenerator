// Draw di-electron invariant mass

void SimpleDraw(const char *inputFileList)
{
  TChain chain("PGS");

  if(!FillChain(chain, inputFileList)) return;

  ExRootResult result;

  TH1 *mass_hist = result.AddHist1D("mass_hist", "di-electron mass", "di-electron mass, GeV/c^{2}", "number of entries", 70, 60.0, 130.0);

  chain.Draw("sqrt((Electron.E[0]+Electron.E[1])^2 - (Electron.Px[0]+Electron.Px[1])^2 - (Electron.Py[0]+Electron.Py[1])^2 - (Electron.Pz[0]+Electron.Pz[1])^2) >> mass_hist", "Electron_size > 1");

  result.Print(); 
}
