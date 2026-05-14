// ============================================================
// reco2018tree.cc
// TTree producer for NA64 2018 calibration data (100 GeV beam)
// Runs of interest: 3610, 3655, 3661, 3662, 3663, 3664, 3665
//
// Key feature: Hod2018 (HOD2) raw SADC waveform access.
// SADC preprocessing mode 0-33 stores the full waveform.
// Standard RunCellReco fails (no 2018 HOD2 calibration loaded).
// Solution: read raw ChipSADC::Digit samples directly and compute
//   amplitude = max(samples) - pedestal(mean of first 5 samples)
//
// Detectors with data for these runs:
//   ECAL, HCAL3, SRD, VETO, HOD2
//   (MM, GEM, Straw, BGO, tracking: empty - B field = 0T)
//
// Prerequisites (see README.md):
//   - ENABLE_TRACKING=1, GENFIT_DIR, TTOOLS_DIR, CATSC_DIR exported
//   - conddb.cc run range fix (bug: runs 3574-3854 had no geometry)
//   - led.cc LED fallback fix (bug: run 3661 has no LED calibration)
//   - maps/2018.xml/HODO.xml run range extended to 3610 (bug: was 3667)
// ============================================================

#include "DaqEventsManager.h"
#include "DaqEvent.h"
#include "ChipSADC.h"
using namespace CS;

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
using namespace std;

#include "p348reco.h"
#include "shower.h"
#include "badburst.h"
#include "timecut.h"

int main(int argc, char *argv[])
{
  if (argc < 2) {
    cerr << "Usage: ./reco2018tree.exe file1 [file2 ...]" << endl;
    return 1;
  }

  DaqEventsManager manager;
  manager.SetMapsDir("../maps");
  for (int i = 1; i < argc; ++i)
    manager.AddDataSource(argv[i]);
  manager.Print();

  TFile* file = new TFile("reco2018tree.root", "RECREATE");
  TTree* tree = new TTree("tree", "NA64 2018 reco tree");

  // ===== Event ID =====
  UInt_t run, spill, spillevent, runevent;
  tree->Branch("run",        &run);
  tree->Branch("spill",      &spill);
  tree->Branch("spillevent", &spillevent);
  tree->Branch("runevent",   &runevent);

  // ===== Trigger flags =====
  bool isPhysics, isCalibration, isOnSpill;
  bool isTriggerPhysics, isTriggerBeam, isTriggerRandom, isTriggerECAL, isTriggerMuonInv;
  tree->Branch("isPhysics",        &isPhysics);
  tree->Branch("isCalibration",    &isCalibration);
  tree->Branch("isOnSpill",        &isOnSpill);
  tree->Branch("isTriggerPhysics", &isTriggerPhysics);
  tree->Branch("isTriggerBeam",    &isTriggerBeam);
  tree->Branch("isTriggerRandom",  &isTriggerRandom);
  tree->Branch("isTriggerECAL",    &isTriggerECAL);
  tree->Branch("isTriggerMuonInv", &isTriggerMuonInv);

  // ===== Timing =====
  double masterTime;
  tree->Branch("masterTime", &masterTime);

  // ===== Scintillators =====
  double S0_amp, S1_amp, S2_amp;
  double S0_t0,  S1_t0,  S2_t0;
  tree->Branch("S0_amp", &S0_amp); tree->Branch("S0_t0", &S0_t0);
  tree->Branch("S1_amp", &S1_amp); tree->Branch("S1_t0", &S1_t0);
  tree->Branch("S2_amp", &S2_amp); tree->Branch("S2_t0", &S2_t0);

  // ===== ECAL =====
  double ecalTotal, ecal0, ecal1;
  double ECAL_ene[2][6][6], ECAL_t0[2][6][6];
  tree->Branch("ecalTotal", &ecalTotal);
  tree->Branch("ecal0",     &ecal0);
  tree->Branch("ecal1",     &ecal1);
  tree->Branch("ECAL_ene",  ECAL_ene, "ECAL_ene[2][6][6]/D");
  tree->Branch("ECAL_t0",   ECAL_t0,  "ECAL_t0[2][6][6]/D");

  // ===== HCAL =====
  // Note: hcal0/1/2 return NaN for calibration runs (B=0T, no calibration)
  // hcal3 has data (zero-degree calorimeter, separate from main magnet)
  double hcal0, hcal1, hcal2, hcal3;
  tree->Branch("hcal0", &hcal0);
  tree->Branch("hcal1", &hcal1);
  tree->Branch("hcal2", &hcal2);
  tree->Branch("hcal3", &hcal3);

  // ===== SRD =====
  double SRD_ene[4], SRD_t0[4], srdTotal;
  tree->Branch("srdTotal", &srdTotal);
  tree->Branch("SRD_ene",  SRD_ene, "SRD_ene[4]/D");
  tree->Branch("SRD_t0",   SRD_t0,  "SRD_t0[4]/D");

  // ===== VETO =====
  double VETO_ene[6], VETO_t0[6];
  tree->Branch("VETO_ene", VETO_ene, "VETO_ene[6]/D");
  tree->Branch("VETO_t0",  VETO_t0,  "VETO_t0[6]/D");

  // ===== ECAL shower chi2 =====
  double showerChi2;
  tree->Branch("showerChi2", &showerChi2);

  // ================================================================
  // HOD2 (Hod2018): 16 X strips + 16 Y strips
  //
  // Two sets of branches:
  //
  // 1) RAW: amplitude computed directly from ChipSADC waveform samples
  //    amp = max(samples) - pedestal
  //    pedestal = mean of first 5 samples
  //    → This is the primary source of HOD2 data for these runs
  //
  // 2) RECO: output of standard RunCellReco (Cell.amplitude, Cell.energy)
  //    → Also has data after fixing the maps run range (HODO.xml)
  //    → Kept for cross-checking with raw amplitudes
  //
  // Channel mapping (HODO.xml, port 3, source 623):
  //   HOD2X: ADC channels 32-47 → strips 15-0 (reversed)
  //   HOD2Y: ADC channels 48-63 → strips 15-0 (reversed)
  // Gains (from Ulloa 2019 report, run 3661):
  //   X strips: ~9-17 mV/pixel
  //   Y strips: ~9-11 mV/pixel
  // ================================================================

  // HOD2 raw (from ChipSADC::GetSamples)
  double HOD2_xamp[16],      HOD2_yamp[16];      // signal = max - pedestal
  double HOD2_xped[16],      HOD2_yped[16];      // pedestal (mean first 5 samples)
  double HOD2_xmax[16],      HOD2_ymax[16];      // raw max sample value
  int    HOD2_xnsamples[16], HOD2_ynsamples[16]; // waveform length
  tree->Branch("HOD2_xamp",      HOD2_xamp,      "HOD2_xamp[16]/D");
  tree->Branch("HOD2_yamp",      HOD2_yamp,      "HOD2_yamp[16]/D");
  tree->Branch("HOD2_xped",      HOD2_xped,      "HOD2_xped[16]/D");
  tree->Branch("HOD2_yped",      HOD2_yped,      "HOD2_yped[16]/D");
  tree->Branch("HOD2_xmax",      HOD2_xmax,      "HOD2_xmax[16]/D");
  tree->Branch("HOD2_ymax",      HOD2_ymax,      "HOD2_ymax[16]/D");
  tree->Branch("HOD2_xnsamples", HOD2_xnsamples, "HOD2_xnsamples[16]/I");
  tree->Branch("HOD2_ynsamples", HOD2_ynsamples, "HOD2_ynsamples[16]/I");

  // HOD2 standard reco (from RunCellReco via RunP348Reco)
  double HOD2_xamp_reco[16], HOD2_yamp_reco[16];
  double HOD2_xene_reco[16], HOD2_yene_reco[16];
  double HOD2_xt0_reco[16],  HOD2_yt0_reco[16];
  bool   HOD2_xhit_reco[16], HOD2_yhit_reco[16];
  vector<int> hod2_xhits_reco, hod2_yhits_reco;
  tree->Branch("HOD2_xamp_reco", HOD2_xamp_reco, "HOD2_xamp_reco[16]/D");
  tree->Branch("HOD2_yamp_reco", HOD2_yamp_reco, "HOD2_yamp_reco[16]/D");
  tree->Branch("HOD2_xene_reco", HOD2_xene_reco, "HOD2_xene_reco[16]/D");
  tree->Branch("HOD2_yene_reco", HOD2_yene_reco, "HOD2_yene_reco[16]/D");
  tree->Branch("HOD2_xt0_reco",  HOD2_xt0_reco,  "HOD2_xt0_reco[16]/D");
  tree->Branch("HOD2_yt0_reco",  HOD2_yt0_reco,  "HOD2_yt0_reco[16]/D");
  tree->Branch("HOD2_xhit_reco", HOD2_xhit_reco, "HOD2_xhit_reco[16]/O");
  tree->Branch("HOD2_yhit_reco", HOD2_yhit_reco, "HOD2_yhit_reco[16]/O");
  tree->Branch("hod2_xhits_reco", &hod2_xhits_reco);
  tree->Branch("hod2_yhits_reco", &hod2_yhits_reco);

  print_software_version_banner();

  // ==================== EVENT LOOP ====================
  while (manager.ReadEvent()) {
    const int nevt = manager.GetEventsCounter();
    if (nevt % 1000 == 1)
      cout << "===> Event #" << nevt << endl;

    if (!manager.DecodeEvent()) {
      cout << "WARNING: fail to decode event #" << nevt << endl;
      continue;
    }

    RecoEvent e = RunP348Reco(manager);
    timecut(e, 4.);

    // --- Event ID ---
    run        = e.run;        spill      = e.spill;
    spillevent = e.spillevent; runevent   = e.runevent;

    // --- Trigger ---
    isPhysics        = e.isPhysics;      isCalibration    = e.isCalibration;
    isOnSpill        = e.isOnSpill;      isTriggerPhysics = e.isTriggerPhysics;
    isTriggerBeam    = e.isTriggerBeam;  isTriggerRandom  = e.isTriggerRandom;
    isTriggerECAL    = e.isTriggerECAL;  isTriggerMuonInv = e.isTriggerMuonInv;

    // --- Timing ---
    masterTime = e.masterTime;

    // --- Scintillators ---
    S0_amp = e.S0.energy; S0_t0 = e.S0.t0ns();
    S1_amp = e.S1.energy; S1_t0 = e.S1.t0ns();
    S2_amp = e.S2.energy; S2_t0 = e.S2.t0ns();

    // --- ECAL ---
    ecal0     = e.ecalTotalEnergy(0);
    ecal1     = e.ecalTotalEnergy(1);
    ecalTotal = ecal0 + ecal1;
    for (int d = 0; d < 2; d++)
      for (int x = 0; x < 6; x++)
        for (int y = 0; y < 6; y++) {
          ECAL_ene[d][x][y] = e.ECAL[d][x][y].energy;
          ECAL_t0[d][x][y]  = e.ECAL[d][x][y].t0ns();
        }

    // --- HCAL ---
    hcal0 = e.hcalTotalEnergy(0); hcal1 = e.hcalTotalEnergy(1);
    hcal2 = e.hcalTotalEnergy(2); hcal3 = e.hcalTotalEnergy(3);

    // --- SRD ---
    srdTotal = 0;
    for (int i = 0; i < 4; i++) {
      SRD_ene[i] = e.SRD[i].energy;
      SRD_t0[i]  = e.SRD[i].t0ns();
      srdTotal  += e.SRD[i].energy;
    }

    // --- VETO ---
    for (int i = 0; i < 6; i++) {
      VETO_ene[i] = e.VETO[i].energy;
      VETO_t0[i]  = e.VETO[i].t0ns();
    }

    // --- Shower chi2 ---
    showerChi2 = calcShowerChi2(e);

    // -------------------------------------------------------
    // HOD2 RAW: read waveform directly from ChipSADC digits
    // Bypasses RunCellReco (fails for preprocessing mode 0-33)
    // -------------------------------------------------------
    for (int i = 0; i < 16; i++) {
      HOD2_xamp[i] = HOD2_yamp[i] = 0;
      HOD2_xped[i] = HOD2_yped[i] = 0;
      HOD2_xmax[i] = HOD2_ymax[i] = 0;
      HOD2_xnsamples[i] = HOD2_ynsamples[i] = 0;
    }

    for (auto it = manager.GetEventDigits().begin();
              it != manager.GetEventDigits().end(); ++it) {
      const CS::ChipSADC::Digit* sadc =
        dynamic_cast<const CS::ChipSADC::Digit*>(it->second);
      if (!sadc) continue;

      const string dname = sadc->GetDetID().GetName();
      if (dname != "HOD2X" && dname != "HOD2Y") continue;

      const int strip = sadc->GetY();
      if (strip < 0 || strip > 15) continue;

      const vector<CS::uint16>& samples = sadc->GetSamples();
      if (samples.empty()) continue;

      const int ns = (int)samples.size();

      // Pedestal = mean of first 5 samples
      int nped = min(ns, 5);
      double ped = 0;
      for (int s = 0; s < nped; s++) ped += samples[s];
      ped /= nped;

      CS::uint16 maxS = *max_element(samples.begin(), samples.end());
      double amp = (double)maxS - ped;

      if (dname == "HOD2X") {
        HOD2_xamp[strip]      = amp;
        HOD2_xped[strip]      = ped;
        HOD2_xmax[strip]      = (double)maxS;
        HOD2_xnsamples[strip] = ns;
      } else {
        HOD2_yamp[strip]      = amp;
        HOD2_yped[strip]      = ped;
        HOD2_ymax[strip]      = (double)maxS;
        HOD2_ynsamples[strip] = ns;
      }
    }

    // -------------------------------------------------------
    // HOD2 RECO: standard RunCellReco output (cross-check)
    // -------------------------------------------------------
    for (int i = 0; i < 16; i++) {
      HOD2_xamp_reco[i] = e.HOD2.xplane[i].amplitude;
      HOD2_xene_reco[i] = e.HOD2.xplane[i].energy;
      HOD2_xt0_reco[i]  = e.HOD2.xplane[i].t0ns();
      HOD2_xhit_reco[i] = e.HOD2.xplane[i].hasDigit;
      HOD2_yamp_reco[i] = e.HOD2.yplane[i].amplitude;
      HOD2_yene_reco[i] = e.HOD2.yplane[i].energy;
      HOD2_yt0_reco[i]  = e.HOD2.yplane[i].t0ns();
      HOD2_yhit_reco[i] = e.HOD2.yplane[i].hasDigit;
    }
    hod2_xhits_reco = e.HOD2.xhits;
    hod2_yhits_reco = e.HOD2.yhits;

    tree->Fill();
  }
  // ==================== END EVENT LOOP ====================

  cout << "Writing output to reco2018tree.root" << endl;
  file->Write();
  file->Close();
  return 0;
}
