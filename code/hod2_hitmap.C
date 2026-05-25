// This macros is useful for plotting amplitude vs channel for the hodoscope. It will also plot the 31x31 distribution in the X-Y plane
// It will assign the position where the biggest amplitude is among vecinity.

{
  double xamp_p3[16], yamp_p3[16];
  tree->SetBranchAddress("HOD2_xamp_p3", xamp_p3);
  tree->SetBranchAddress("HOD2_yamp_p3", yamp_p3);

  TH2F* h2 = new TH2F("h2","HOD2 hit map;X (mm);Y (mm)",
                        33, 0, 33, 33, 0, 33);
  TH2F* hXamp = new TH2F("hXamp","HOD2X amp vs channel;channel;amplitude",
                           16, 0, 16, 100, 0, 500);
  TH2F* hYamp = new TH2F("hYamp","HOD2Y amp vs channel;channel;amplitude",
                           16, 0, 16, 100, 0, 500);

  double threshold = 50;

  for(Long64_t ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);

    int x_strip = -1; double x_max = threshold;
    for(int i=0; i<16; i++)
      if(xamp_p3[i] > x_max){ x_max = xamp_p3[i]; x_strip = i; }

    double x_mm = -1;
    if(x_strip >= 0){
      bool right = (x_strip < 15) && (xamp_p3[x_strip+1] > threshold);
      bool left  = (x_strip >  0) && (xamp_p3[x_strip-1] > threshold);
      if(right && !left) x_mm = 2*(x_strip+1);  
      else if(left && !right) x_mm = 2*x_strip;  
      else x_mm = 2*x_strip + 1;                
      x_mm = max(1.0, min(31.0, x_mm));
    }

    int y_strip = -1; double y_max = threshold;
    for(int i=0; i<16; i++)
      if(yamp_p3[i] > y_max){ y_max = yamp_p3[i]; y_strip = i; }

    double y_mm = -1;
    if(y_strip >= 0){
      bool right = (y_strip < 15) && (yamp_p3[y_strip+1] > threshold);
      bool left  = (y_strip >  0) && (yamp_p3[y_strip-1] > threshold);
      if(right && !left) y_mm = 2*(y_strip+1);
      else if(left && !right) y_mm = 2*y_strip;
      else y_mm = 2*y_strip + 1;
      y_mm = max(1.0, min(31.0, y_mm));
    }

    if(x_mm > 0 && y_mm > 0)
      h2->Fill(x_mm, y_mm);

    // Amp vs canal (port 3: canal = 15-strip)
    for(int i=0; i<16; i++){
      if(xamp_p3[i] > 0) hXamp->Fill(15-i, xamp_p3[i]);
      if(yamp_p3[i] > 0) hYamp->Fill(15-i, yamp_p3[i]);
    }
  }

  TCanvas* c1 = new TCanvas("c1","Hit map",800,600);
  gStyle->SetOptStat(1);
  h2->Draw("COLZ");

  TCanvas* c2 = new TCanvas("c2","Amp vs channel",1200,500);
  c2->Divide(2,1);
  c2->cd(1); hXamp->Draw("COL");
  c2->cd(2); hYamp->Draw("COL");
}
