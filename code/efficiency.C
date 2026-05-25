{
  double xamp_p3[16], yamp_p3[16];
  double S2_amp;
  tree->SetBranchAddress("HOD2_xamp_p3", xamp_p3);
  tree->SetBranchAddress("HOD2_yamp_p3", yamp_p3);
  tree->SetBranchAddress("S2_amp", &S2_amp);

  double threshold_S2  = 700;  // ADC units (calibrated)
  double threshold_HOD = 50;   // ADC units (raw waveform)

  Long64_t n_S2       = 0;  // events with S2
  Long64_t n_S2_HOD2  = 0;  // eventos with S2 and HOD2

  for(Long64_t ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);

    // S2 trigger
    if(S2_amp < threshold_S2) continue;
    n_S2++;

    // gives AND condition
    bool hod2_x_hit = false, hod2_y_hit = false;
    for(int i=0; i<16; i++){
      if(xamp_p3[i] > threshold_HOD) hod2_x_hit = true;
      if(yamp_p3[i] > threshold_HOD) hod2_y_hit = true;
    }

    if(hod2_x_hit && hod2_y_hit) n_S2_HOD2++;
  }

  double efficiency = (n_S2 > 0) ? (double)n_S2_HOD2 / n_S2 : 0;

  cout << "===== HOD2 Efficiency =====" << endl;
  cout << "Events with S2 (beam):     " << n_S2      << endl;
  cout << "Events with S2 AND HOD2:   " << n_S2_HOD2 << endl;
  cout << "Eficiency:                " << efficiency*100 << " %" << endl;
  cout << "===========================" << endl;
}
