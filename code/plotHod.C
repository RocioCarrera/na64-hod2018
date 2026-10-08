void plotHod() {

  TFile* f = new TFile("reco2018tree.root");
  TTree* tree = (TTree*)f->Get("tree");

  // Events with HOD2 signal
  cout << "Events with HOD2_xamp_p3[0]>0: " << tree->GetEntries("HOD2_xamp_p3[0]>0") << endl;

  // Amplitude distribution for strip X0
  new TCanvas("c0","Strip X0");
  tree->Draw("HOD2_xamp_p3[0]", "HOD2_xamp_p3[0]>5");

  // Beam profile: hits per X strip
  TH1F* hx = new TH1F("hx","Beam profile X;Strip;Events",16,0,16);
  for(int i=0; i<16; i++)
    hx->SetBinContent(i+1, tree->GetEntries(Form("HOD2_xamp_p3[%d]>50",i)));
  new TCanvas("c1","Beam profile X");
  hx->SetFillColor(kAzure-9);
  hx->Draw();

  // 2D hit map (strip X vs strip Y)
  TH2F* h2 = new TH2F("h2","HOD2 hit map;X strip;Y strip",16,0,16,16,0,16);
  double xamp[16], yamp[16];
  tree->SetBranchAddress("HOD2_xamp_p3", xamp);
  tree->SetBranchAddress("HOD2_yamp_p3", yamp);
  for(Long64_t ev=0; ev<tree->GetEntries(); ev++) {
    tree->GetEntry(ev);
    int xstrip=-1; double xmax=50;
    for(int i=0;i<16;i++) if(xamp[i]>xmax){xmax=xamp[i]; xstrip=i;}
    int ystrip=-1; double ymax=50;
    for(int i=0;i<16;i++) if(yamp[i]>ymax){ymax=yamp[i]; ystrip=i;}
    if(xstrip>=0 && ystrip>=0) h2->Fill(xstrip,ystrip);
  }
  new TCanvas("c2","Hit map");
  h2->SetDirectory(0);
  h2->Draw("COLZ");
}
