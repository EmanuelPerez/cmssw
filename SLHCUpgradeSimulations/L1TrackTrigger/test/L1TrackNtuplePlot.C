// ----------------------------------------------------------------------------------------------------------------
// Basic example script for making tracking performance plots using the ntuples produced by L1TrackNtupleMaker.cc
// By Louise Skinnari, June 2013  
// ----------------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text); 


// ----------------------------------------------------------------------------------------------------------------
// Main script
void L1TrackNtuplePlot(TString type) {

  SetPlotStyle();
  gROOT->ProcessLine("#include <vector>");

  bool doDR = false; //turn on to use dR matching instead of simtrack ID matching for electrons/pions
  bool doDetailedPlots = false; //turn on to make full set of plots


  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  if (type=="SingleMuon") {
    tree->Add("RootFiles/SingleMuPlus_BE5D_TrkPerf.root");
    tree->Add("RootFiles/SingleMuMinus_BE5D_TrkPerf.root");
  }
  else if (type=="SingleElectron") {
    tree->Add("RootFiles/SingleElectron_BE5D_TrkPerf.root");
    tree->Add("RootFiles/SinglePositron_BE5D_TrkPerf.root");
  }
  else {
    tree->Add("RootFiles/"+type+"_BE5D_TrkPerf.root");
  }

  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }


  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches

  // sim track properties
  vector<float>* simtrk_pt;
  vector<float>* simtrk_eta;
  vector<float>* simtrk_phi;
  vector<float>* simtrk_z0;
  vector<int>*   simtrk_id;   //simtrackID
  vector<int>*   simtrk_type; //PDG ID of the sim track (13==muon, 11==electron, and so on)

  // *sim track* properties, for sim tracks that are matched to an L1 track using simtrackID
  vector<float>* matchID_simtrk_pt;
  vector<float>* matchID_simtrk_eta;
  vector<float>* matchID_simtrk_phi;
  vector<float>* matchID_simtrk_z0;
  vector<int>*   matchID_simtrk_id;
  vector<int>*   matchID_simtrk_type;

  // *L1 track* properties, for sim tracks that are matched to an L1 track using simtrackID
  vector<float>* matchID_trk_pt;
  vector<float>* matchID_trk_eta;
  vector<float>* matchID_trk_phi;
  vector<float>* matchID_trk_z0;
  vector<float>* matchID_trk_chi2; 
  vector<int>*   matchID_trk_nstub;
  vector<int>*   matchID_trk_nmatch;

  // *sim track* properties, for sim tracks that are matched to an L1 track using dR
  vector<float>* matchDR_simtrk_pt;
  vector<float>* matchDR_simtrk_eta;
  vector<float>* matchDR_simtrk_phi;
  vector<float>* matchDR_simtrk_z0;
  vector<int>*   matchDR_simtrk_id;
  vector<int>*   matchDR_simtrk_type;

  // *L1 track* properties, for sim tracks that are matched to an L1 track using dR
  vector<float>* matchDR_trk_pt;
  vector<float>* matchDR_trk_eta;
  vector<float>* matchDR_trk_phi;
  vector<float>* matchDR_trk_z0;
  vector<float>* matchDR_trk_chi2; 
  vector<int>*   matchDR_trk_nstub;
  vector<int>*   matchDR_trk_nmatch;

  
  TBranch* b_simtrk_pt;
  TBranch* b_simtrk_eta;
  TBranch* b_simtrk_phi;
  TBranch* b_simtrk_z0;
  TBranch* b_simtrk_id;
  TBranch* b_simtrk_type;

  TBranch* b_matchID_simtrk_pt;
  TBranch* b_matchID_simtrk_eta;
  TBranch* b_matchID_simtrk_phi;
  TBranch* b_matchID_simtrk_z0;
  TBranch* b_matchID_simtrk_id;
  TBranch* b_matchID_simtrk_type;

  TBranch* b_matchID_trk_pt;
  TBranch* b_matchID_trk_eta;
  TBranch* b_matchID_trk_phi;
  TBranch* b_matchID_trk_z0;
  TBranch* b_matchID_trk_chi2; 
  TBranch* b_matchID_trk_nstub;
  TBranch* b_matchID_trk_nmatch;

  TBranch* b_matchDR_simtrk_pt;
  TBranch* b_matchDR_simtrk_eta;
  TBranch* b_matchDR_simtrk_phi;
  TBranch* b_matchDR_simtrk_z0;
  TBranch* b_matchDR_simtrk_id;
  TBranch* b_matchDR_simtrk_type;

  TBranch* b_matchDR_trk_pt;
  TBranch* b_matchDR_trk_eta;
  TBranch* b_matchDR_trk_phi;
  TBranch* b_matchDR_trk_z0;
  TBranch* b_matchDR_trk_chi2; 
  TBranch* b_matchDR_trk_nstub;
  TBranch* b_matchDR_trk_nmatch;


  simtrk_pt   = 0;
  simtrk_eta  = 0;
  simtrk_phi  = 0;
  simtrk_z0   = 0;
  simtrk_id   = 0;
  simtrk_type = 0;

  matchID_simtrk_pt  = 0;
  matchID_simtrk_eta = 0;
  matchID_simtrk_phi = 0;
  matchID_simtrk_z0  = 0;
  matchID_simtrk_id  = 0;
  matchID_simtrk_type  = 0;

  matchID_trk_pt  = 0;
  matchID_trk_eta = 0;
  matchID_trk_phi = 0;
  matchID_trk_z0  = 0;
  matchID_trk_chi2   = 0; 
  matchID_trk_nstub  = 0;
  matchID_trk_nmatch = 0;
  
  matchDR_simtrk_pt  = 0;
  matchDR_simtrk_eta = 0;
  matchDR_simtrk_phi = 0;
  matchDR_simtrk_z0  = 0;
  matchDR_simtrk_id  = 0;
  matchDR_simtrk_type  = 0;

  matchDR_trk_pt  = 0;
  matchDR_trk_eta = 0;
  matchDR_trk_phi = 0;
  matchDR_trk_z0  = 0;
  matchDR_trk_chi2   = 0; 
  matchDR_trk_nstub  = 0;
  matchDR_trk_nmatch = 0;


  tree->SetBranchAddress("simtrk_pt",  &simtrk_pt,  &b_simtrk_pt);
  tree->SetBranchAddress("simtrk_eta", &simtrk_eta, &b_simtrk_eta);
  tree->SetBranchAddress("simtrk_phi", &simtrk_phi, &b_simtrk_phi);
  tree->SetBranchAddress("simtrk_z0",  &simtrk_z0,  &b_simtrk_z0);
  tree->SetBranchAddress("simtrk_id",  &simtrk_id,  &b_simtrk_id);
  tree->SetBranchAddress("simtrk_type",&simtrk_type,&b_simtrk_type);

  tree->SetBranchAddress("matchID_simtrk_pt",  &matchID_simtrk_pt,  &b_matchID_simtrk_pt);
  tree->SetBranchAddress("matchID_simtrk_eta", &matchID_simtrk_eta, &b_matchID_simtrk_eta);
  tree->SetBranchAddress("matchID_simtrk_phi", &matchID_simtrk_phi, &b_matchID_simtrk_phi);
  tree->SetBranchAddress("matchID_simtrk_z0",  &matchID_simtrk_z0,  &b_matchID_simtrk_z0);
  tree->SetBranchAddress("matchID_simtrk_id",  &matchID_simtrk_id,  &b_matchID_simtrk_id);
  tree->SetBranchAddress("matchID_simtrk_type",&matchID_simtrk_type,&b_matchID_simtrk_type);

  tree->SetBranchAddress("matchID_trk_pt",    &matchID_trk_pt,    &b_matchID_trk_pt);
  tree->SetBranchAddress("matchID_trk_eta",   &matchID_trk_eta,   &b_matchID_trk_eta);
  tree->SetBranchAddress("matchID_trk_phi",   &matchID_trk_phi,   &b_matchID_trk_phi);
  tree->SetBranchAddress("matchID_trk_z0",    &matchID_trk_z0,    &b_matchID_trk_z0);
  tree->SetBranchAddress("matchID_trk_chi2",  &matchID_trk_chi2,  &b_matchID_trk_chi2);
  tree->SetBranchAddress("matchID_trk_nstub", &matchID_trk_nstub, &b_matchID_trk_nstub);
  tree->SetBranchAddress("matchID_trk_nmatch",&matchID_trk_nmatch,&b_matchID_trk_nmatch);

  tree->SetBranchAddress("matchDR_simtrk_pt",  &matchDR_simtrk_pt,  &b_matchDR_simtrk_pt);
  tree->SetBranchAddress("matchDR_simtrk_eta", &matchDR_simtrk_eta, &b_matchDR_simtrk_eta);
  tree->SetBranchAddress("matchDR_simtrk_phi", &matchDR_simtrk_phi, &b_matchDR_simtrk_phi);
  tree->SetBranchAddress("matchDR_simtrk_z0",  &matchDR_simtrk_z0,  &b_matchDR_simtrk_z0);
  tree->SetBranchAddress("matchDR_simtrk_id",  &matchDR_simtrk_id,  &b_matchDR_simtrk_id);
  tree->SetBranchAddress("matchDR_simtrk_type",&matchDR_simtrk_type,&b_matchDR_simtrk_type);

  tree->SetBranchAddress("matchDR_trk_pt",    &matchDR_trk_pt,    &b_matchDR_trk_pt);
  tree->SetBranchAddress("matchDR_trk_eta",   &matchDR_trk_eta,   &b_matchDR_trk_eta);
  tree->SetBranchAddress("matchDR_trk_phi",   &matchDR_trk_phi,   &b_matchDR_trk_phi);
  tree->SetBranchAddress("matchDR_trk_z0",    &matchDR_trk_z0,    &b_matchDR_trk_z0);
  tree->SetBranchAddress("matchDR_trk_chi2",  &matchDR_trk_chi2,  &b_matchDR_trk_chi2);
  tree->SetBranchAddress("matchDR_trk_nstub", &matchDR_trk_nstub, &b_matchDR_trk_nstub);
  tree->SetBranchAddress("matchDR_trk_nmatch",&matchDR_trk_nmatch,&b_matchDR_trk_nmatch);
  

  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  /////////////////////////////////////////////////
  // NOTATION:                                   //
  // 'C' - Central eta range, |eta|<0.8          //
  // 'I' - Intermediate eta range, 0.8<|eta|<1.6 //
  // 'F' - Forward eta range, |eta|>1.6          //
  //                                             //
  // 'L' - Low pt range, pt<5 GeV                //
  // 'M' - Middle pt range, 5<pt<15 GeV          //
  // 'H' - High pt range, pt>15 GeV              //
  /////////////////////////////////////////////////


  TH1F* h_simtrk_pt   = new TH1F("simtrk_pt",   ";Sim track p_{T} [GeV]; Sim tracks / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_simtrk_pt_L = new TH1F("simtrk_pt_L", ";Sim track p_{T} [GeV]; Sim tracks / 0.1 GeV",  50,  0,     5.0);
  TH1F* h_simtrk_eta  = new TH1F("simtrk_eta",  ";Sim track #eta; Sim tracks / 0.1",             50, -2.5,   2.5);
  TH1F* h_simtrk_phi  = new TH1F("simtrk_phi",  ";Sim track #phi [rad]; Sim tracks / 0.1",       64, -3.2,   3.2);

  TH1F* h_matchID_simtrk_pt   = new TH1F("matchID_simtrk_pt",   ";Sim track p_{T} [GeV]; Sim tracks / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_matchID_simtrk_pt_L = new TH1F("matchID_simtrk_pt_L", ";Sim track p_{T} [GeV]; Sim tracks / 0.1 GeV",  50,  0,     5.0);
  TH1F* h_matchID_simtrk_eta  = new TH1F("matchID_simtrk_eta",  ";Sim track #eta; Sim tracks / 0.1",             50, -2.5,   2.5);
  TH1F* h_matchID_simtrk_phi  = new TH1F("matchID_simtrk_phi",  ";Sim track #phi [rad]; Sim tracks / 0.1",       64, -3.2,   3.2);

  TH1F* h_matchID_trk_nstub   = new TH1F("matchID_trk_nstub",   ";Number of stubs; L1 tracks / 1.0", 15, 0, 15);
  TH1F* h_matchID_trk_nstub_C = new TH1F("matchID_trk_nstub_C", ";Number of stubs; L1 tracks / 1.0", 15, 0, 15);
  TH1F* h_matchID_trk_nstub_I = new TH1F("matchID_trk_nstub_I", ";Number of stubs; L1 tracks / 1.0", 15, 0, 15);
  TH1F* h_matchID_trk_nstub_F = new TH1F("matchID_trk_nstub_F", ";Number of stubs; L1 tracks / 1.0", 15, 0, 15);

  // chi2 histograms
  // note: last bin is an overflow bin
  TH1F* h_matchID_trk_chi2     = new TH1F("matchID_trk_chi2",     ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_C_L = new TH1F("matchID_trk_chi2_C_L", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_I_L = new TH1F("matchID_trk_chi2_I_L", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_F_L = new TH1F("matchID_trk_chi2_F_L", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_C_M = new TH1F("matchID_trk_chi2_C_M", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_I_M = new TH1F("matchID_trk_chi2_I_M", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_F_M = new TH1F("matchID_trk_chi2_F_M", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_C_H = new TH1F("matchID_trk_chi2_C_H", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_I_H = new TH1F("matchID_trk_chi2_I_H", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);
  TH1F* h_matchID_trk_chi2_F_H = new TH1F("matchID_trk_chi2_F_H", ";#chi^{2}; L1 tracks / 1.0", 100, 0, 100);

  // chi2/dof histograms
  // note: lastbin is an overflow bin
  TH1F* h_matchID_trk_chi2_dof     = new TH1F("matchID_trk_chi2_dof",     ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_C_L = new TH1F("matchID_trk_chi2_dof_C_L", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_I_L = new TH1F("matchID_trk_chi2_dof_I_L", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_F_L = new TH1F("matchID_trk_chi2_dof_F_L", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_C_M = new TH1F("matchID_trk_chi2_dof_C_M", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_I_M = new TH1F("matchID_trk_chi2_dof_I_M", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_F_M = new TH1F("matchID_trk_chi2_dof_F_M", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_C_H = new TH1F("matchID_trk_chi2_dof_C_H", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_I_H = new TH1F("matchID_trk_chi2_dof_I_H", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);
  TH1F* h_matchID_trk_chi2_dof_F_H = new TH1F("matchID_trk_chi2_dof_F_H", ";#chi^{2} / D.O.F.; L1 tracks / 0.1", 150, 0, 15);

  // resolution histograms
  TH1F* h_res_pt    = new TH1F("res_pt",    ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.05",   200,-5.0,   5.0);
  TH1F* h_res_ptRel = new TH1F("res_ptRel", ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.01", 200,-1.0,   1.0);
  TH1F* h_res_eta   = new TH1F("res_eta",   ";#eta residual (L1 - sim); L1 tracks / 0.0002",        100,-0.01, 0.01);
  TH1F* h_res_phi   = new TH1F("res_phi",   ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001",  100,-0.005,0.005);
  TH1F* h_res_z0    = new TH1F("res_z0",    ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02",    100,-1,     1);
  TH1F* h_res_z0_C  = new TH1F("res_z0_C",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02",    100,-1,     1);
  TH1F* h_res_z0_I  = new TH1F("res_z0_I",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02",    100,-1,     1);
  TH1F* h_res_z0_F  = new TH1F("res_z0_F",  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02",    100,-1,     1);

  // resolution vs. pt histograms
  const int nRANGE = 20;
  TString ptrange[nRANGE] = {"0-5","5-10", "10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50-55",
  			     "55-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"};

  TH1F* h_resVsPt_pt[nRANGE];
  TH1F* h_resVsPt_pt_C[nRANGE];
  TH1F* h_resVsPt_pt_I[nRANGE];
  TH1F* h_resVsPt_pt_F[nRANGE];

  TH1F* h_resVsPt_ptRel[nRANGE];
  TH1F* h_resVsPt_ptRel_C[nRANGE];
  TH1F* h_resVsPt_ptRel_I[nRANGE];
  TH1F* h_resVsPt_ptRel_F[nRANGE];

  TH1F* h_resVsPt_z0[nRANGE];
  TH1F* h_resVsPt_z0_C[nRANGE];
  TH1F* h_resVsPt_z0_I[nRANGE];
  TH1F* h_resVsPt_z0_F[nRANGE];

  TH1F* h_resVsPt_phi[nRANGE];
  TH1F* h_resVsPt_phi_C[nRANGE];
  TH1F* h_resVsPt_phi_I[nRANGE];
  TH1F* h_resVsPt_phi_F[nRANGE];

  TH1F* h_resVsPt_eta[nRANGE];

  for (int i=0; i<nRANGE; i++) {
    h_resVsPt_pt[i]   = new TH1F("resVsPt_pt_"+ptrange[i],   ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);
    h_resVsPt_pt_C[i] = new TH1F("resVsPt_pt_C_"+ptrange[i], ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);
    h_resVsPt_pt_I[i] = new TH1F("resVsPt_pt_I_"+ptrange[i], ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);
    h_resVsPt_pt_F[i] = new TH1F("resVsPt_pt_F_"+ptrange[i], ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1", 100, -5.0, 5.0);

    // restictive range: -0.15 to 0.15
    h_resVsPt_ptRel[i]   = new TH1F("resVsPt_ptRel_"+ptrange[i],   ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);
    h_resVsPt_ptRel_C[i] = new TH1F("resVsPt_ptRel_c_"+ptrange[i], ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);
    h_resVsPt_ptRel_I[i] = new TH1F("resVsPt_ptRel_I_"+ptrange[i], ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);
    h_resVsPt_ptRel_F[i] = new TH1F("resVsPt_ptRel_F_"+ptrange[i], ";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 300, -0.15, 0.15);

    h_resVsPt_z0[i]   = new TH1F("resVsPt_z0_"+ptrange[i],   ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);
    h_resVsPt_z0_C[i] = new TH1F("resVsPt_z0_C_"+ptrange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);
    h_resVsPt_z0_I[i] = new TH1F("resVsPt_z0_I_"+ptrange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);
    h_resVsPt_z0_F[i] = new TH1F("resVsPt_z0_F_"+ptrange[i], ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.02", 100, -1, 1);

    h_resVsPt_phi[i]   = new TH1F("resVsPt_phi_"+ptrange[i],   ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);
    h_resVsPt_phi_C[i] = new TH1F("resVsPt_phi_C_"+ptrange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);
    h_resVsPt_phi_I[i] = new TH1F("resVsPt_phi_I_"+ptrange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);
    h_resVsPt_phi_F[i] = new TH1F("resVsPt_phi_F_"+ptrange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100, -0.005, 0.005);

    h_resVsPt_eta[i] = new TH1F("resVsPt_eta_"+ptrange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100, -0.01, 0.01);
  }

  // resolution vs. eta histograms
  const int nETARANGE = 25;
  TString etarange[nETARANGE] = {"0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0",
				 "1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0",
				 "2.1","2.2","2.3","2.4","2.5"};

  TH1F* h_resVsEta_eta[nETARANGE];
  TH1F* h_resVsEta_eta_L[nETARANGE];
  TH1F* h_resVsEta_eta_M[nETARANGE];
  TH1F* h_resVsEta_eta_H[nETARANGE];

  TH1F* h_resVsEta_z0[nETARANGE];
  TH1F* h_resVsEta_phi[nETARANGE];
  TH1F* h_resVsEta_pt[nETARANGE];
  TH1F* h_resVsEta_ptRel[nETARANGE];

  for (int i=0; i<nETARANGE; i++) {
    h_resVsEta_eta[i]   = new TH1F("resVsEta_eta_"+etarange[i],   ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);
    h_resVsEta_eta_L[i] = new TH1F("resVsEta_eta_L_"+etarange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);
    h_resVsEta_eta_M[i] = new TH1F("resVsEta_eta_M_"+etarange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);
    h_resVsEta_eta_H[i] = new TH1F("resVsEta_eta_H_"+etarange[i], ";#eta residual (L1 - sim); L1 tracks / 0.0002", 100,-0.01,0.01);

    h_resVsEta_z0[i]  = new TH1F("resVsEta_z0_"+etarange[i],  ";z_{0} residual (L1 - sim) [cm]; L1 tracks / 0.01",   100,-1,    1);
    h_resVsEta_phi[i] = new TH1F("resVsEta_phi_"+etarange[i], ";#phi residual (L1 - sim) [rad]; L1 tracks / 0.0001", 100,-0.005,0.005);

    h_resVsEta_pt[i]    = new TH1F("resVsEta_pt_"+etarange[i],   ";p_{T} residual (L1 - sim) [GeV]; L1 tracks / 0.1",    100,-5.0,5.0);
    h_resVsEta_ptRel[i] = new TH1F("resVsEta_ptRel_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", 100,-1.0,1.0);
  }

  // 2D histograms
  TH2F* h_2d_logchi2_eta     = new TH2F("2d_logchi2_eta",     ";Sim track #eta; log(#chi^{2})",                    50,-2.5,2.5, 100,-4.0,8.0);
  TH2F* h_2d_logchi2_dof_eta = new TH2F("2d_logchi2_cof_eta", ";Sim track #eta; log(#chi^{2} / D.O.F.)",           50,-2.5,2.5, 100,-4.0,8.0);
  TH2F* h_2d_dz0_eta         = new TH2F("2d_dz0_eta",         ";Sim track #eta; |#Deltaz_{0}| (L1 - sim) [cm]",    50,-2.5,2.5, 100, 0,  1.2);
  TH2F* h_2d_deta_eta        = new TH2F("2d_deta_eta",        ";Sim track #eta; |#Delta#eta| (L1 - sim)",          50,-2.5,2.5, 100, 0,  0.012);
  TH2F* h_2d_dphi_eta        = new TH2F("2d_dphi_eta",        ";Sim track #eta; |#Delta#phi| (L1 - sim) [rad]",    50,-2.5,2.5, 100, 0,  0.007);
  TH2F* h_2d_dptRel_eta      = new TH2F("2d_dptRel_eta",      ";Sim track #eta; |#Deltap_{T} / p_{T}| (L1 - sim)", 50,-2.5,2.5, 100, 0,  1.2);


  
  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  for (int i=0; i<nevt; i++) {

    tree->GetEntry(i,0);
  
    bool foundSim = false;  
    int this_simtrk = -999;
    
    
    // ----------------------------------------------------------------------------------------------------------------
    // sim track loop
    for (int it=0; it<(int)simtrk_pt->size(); it++) {
      
      if (simtrk_pt->at(it) < 0.2) continue;
      if (fabs(simtrk_eta->at(it)) > 2.5) continue;
      if (fabs(simtrk_z0->at(it)) > 30.0) continue;
      
      if (foundSim) {
	cout << "More than one good sim track in event... Ignore and continue..." << endl;
	continue;
      }
      
      foundSim = true;
      this_simtrk = it;
      
      // fill simtrk properties
      h_simtrk_pt->Fill(simtrk_pt->at(it));
      if (simtrk_pt->at(it) < 5.0) h_simtrk_pt_L->Fill(simtrk_pt->at(it));
      
      if (simtrk_pt->at(it) > 2.0) {
	h_simtrk_eta->Fill(simtrk_eta->at(it));
	h_simtrk_phi->Fill(simtrk_phi->at(it));
      }
      
    }//end simtrk loop
    
    
    if (!foundSim) continue;
    
    
    if (doDR && (type.Contains("SingleEl") || type.Contains("SinglePi")) ) {
      matchID_simtrk_pt   = matchDR_simtrk_pt;
      matchID_simtrk_eta  = matchDR_simtrk_eta;
      matchID_simtrk_phi  = matchDR_simtrk_phi;
      matchID_simtrk_z0   = matchDR_simtrk_z0;
      matchID_simtrk_id   = matchDR_simtrk_id;
      matchID_simtrk_type = matchDR_simtrk_type;
      
      matchID_trk_pt     = matchDR_trk_pt;
      matchID_trk_eta    = matchDR_trk_eta;
      matchID_trk_phi    = matchDR_trk_phi;
      matchID_trk_z0     = matchDR_trk_z0;
      matchID_trk_chi2   = matchDR_trk_chi2;
      matchID_trk_nstub  = matchDR_trk_nstub;
      matchID_trk_nmatch = matchDR_trk_nmatch;
    }
    
    
    // ----------------------------------------------------------------------------------------------------------------
    // loop over matched L1tracks-simtracks
    for (int it=0; it<(int)matchID_simtrk_pt->size(); it++) {
      
      if (matchID_simtrk_pt->at(it) < 0.2) continue;
      if (fabs(matchID_simtrk_eta->at(it)) > 2.5) continue;
      if (fabs(matchID_simtrk_z0->at(it)) > 30.0) continue;
      
      // use only tracks with min 4 stubs
      if (matchID_trk_nstub->at(it) < 4) continue;
      
    
      if (it != this_simtrk) {
	cout << "Matched but to wrong sim track, ignore event! (matched to " << it << " instead of " << this_simtrk << ")" << endl;
	continue;
      }
      
      
      // ----------------------------------------------------------------------------------------------------------------
      // fill chi2 & chi2/dof histograms before making chi2 cut
      h_2d_logchi2_eta    ->Fill(matchID_simtrk_eta->at(it), log(matchID_trk_chi2->at(it)));
      h_2d_logchi2_dof_eta->Fill(matchID_simtrk_eta->at(it), log(matchID_trk_chi2->at(it)/(2*matchID_trk_nstub->at(it)-4)));
    
      float chi2 = matchID_trk_chi2->at(it);
      int ndof = 2*matchID_trk_nstub->at(it)-4;
      float chi2dof = (float)chi2/ndof;
      if (chi2 > 100) chi2 = 99.9; //for overflow bin
      if (chi2dof > 15) chi2dof = 14.99; //for overflow bin
  
      h_matchID_trk_chi2->Fill(chi2);
      h_matchID_trk_chi2_dof->Fill(chi2dof);
      
      // central eta
      if (fabs(matchID_trk_eta->at(it)) < 0.8) {
	if (matchID_trk_pt->at(it) < 5) {
	  h_matchID_trk_chi2_C_L->Fill(chi2);
          h_matchID_trk_chi2_dof_C_L->Fill(chi2dof);
	} 
	else if (matchID_trk_pt->at(it) < 15 && matchID_trk_pt->at(it) >= 5) {
	  h_matchID_trk_chi2_C_M->Fill(chi2);
	  h_matchID_trk_chi2_dof_C_M->Fill(chi2dof);
	}
	else {
	  h_matchID_trk_chi2_C_H->Fill(chi2);
	  h_matchID_trk_chi2_dof_C_H->Fill(chi2dof);
	}
      }
      // intermediate eta
      else if (fabs(matchID_trk_eta->at(it)) < 1.6 && fabs(matchID_trk_eta->at(it)) >= 0.8) {
	if (matchID_trk_pt->at(it) < 5) {
          h_matchID_trk_chi2_I_L->Fill(chi2);
          h_matchID_trk_chi2_dof_I_L->Fill(chi2dof);
	}
	else if (matchID_trk_pt->at(it) < 15 && matchID_trk_pt->at(it) >= 5) {
          h_matchID_trk_chi2_I_M->Fill(chi2);
          h_matchID_trk_chi2_dof_I_M->Fill(chi2dof);
	} 
	else {
          h_matchID_trk_chi2_I_H->Fill(chi2);
          h_matchID_trk_chi2_dof_I_H->Fill(chi2dof);
	}
      }
      // forward eta
      else if (fabs(matchID_trk_eta->at(it)) >= 1.6) {
        if (matchID_trk_pt->at(it) < 5) {
          h_matchID_trk_chi2_F_L->Fill(chi2);
          h_matchID_trk_chi2_dof_F_L->Fill(chi2dof);
	} 
	else if (matchID_trk_pt->at(it) < 15 && matchID_trk_pt->at(it) >= 5) {
          h_matchID_trk_chi2_F_M->Fill(chi2);
          h_matchID_trk_chi2_dof_F_M->Fill(chi2dof);
       	} 
	else {
          h_matchID_trk_chi2_F_H->Fill(chi2);
          h_matchID_trk_chi2_dof_F_H->Fill(chi2dof);
	}
      }
      

      // ----------------------------------------------------------------------------------------------------------------
      //cut on chi2
      if (matchID_trk_chi2->at(it) > 100.0) continue;

      
      // ----------------------------------------------------------------------------------------------------------------
      // more plots

      // fill matched track histograms
      h_matchID_simtrk_pt->Fill(matchID_simtrk_pt->at(it));
      if (matchID_simtrk_pt->at(it) < 5) h_matchID_simtrk_pt_L->Fill(matchID_simtrk_pt->at(it));

      if (matchID_simtrk_pt->at(it) > 2) {
	h_matchID_simtrk_eta->Fill(matchID_simtrk_eta->at(it));
	h_matchID_simtrk_phi->Fill(matchID_simtrk_phi->at(it));
      }
      
      // fill nstub histograms
      h_matchID_trk_nstub->Fill(matchID_trk_nstub->at(it));
      if (fabs(matchID_trk_eta->at(it)) < 0.8) h_matchID_trk_nstub_C->Fill(matchID_trk_nstub->at(it));
      else if (fabs(matchID_trk_eta->at(it)) < 1.6 && fabs(matchID_trk_eta->at(it)) >= 0.8) h_matchID_trk_nstub_I->Fill(matchID_trk_nstub->at(it));
      else if (fabs(matchID_trk_eta->at(it)) >= 1.6) h_matchID_trk_nstub_F->Fill(matchID_trk_nstub->at(it));
      

      // fill resolution histograms
      h_res_pt   ->Fill(matchID_trk_pt->at(it)  - matchID_simtrk_pt->at(it));
      h_res_ptRel->Fill((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it));
      h_res_eta  ->Fill(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it));
      h_res_phi  ->Fill(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it));
      h_res_z0   ->Fill(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it));

      
      if (fabs(matchID_trk_eta->at(it)) < 0.8) h_res_z0_C->Fill(matchID_trk_z0->at(it) - matchID_simtrk_z0->at(it));
      else if (fabs(matchID_trk_eta->at(it)) < 1.6 && fabs(matchID_trk_eta->at(it)) >= 0.8) h_res_z0_I->Fill(matchID_trk_z0->at(it) - matchID_simtrk_z0->at(it));
      else if (fabs(matchID_trk_eta->at(it)) >= 1.6) h_res_z0_F->Fill(matchID_trk_z0->at(it) - matchID_simtrk_z0->at(it));
      

      // fill resolution vs. pt histograms    
      for (int im=0; im<nRANGE; im++) {
       if ( (matchID_simtrk_pt->at(it) > (float)im*5.0) && (matchID_simtrk_pt->at(it) < (float)im*5.0+5.0) ) {
	  h_resVsPt_pt[im]   ->Fill(matchID_trk_pt->at(it)  - matchID_simtrk_pt->at(it));
	  h_resVsPt_ptRel[im]->Fill((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it));
	  h_resVsPt_eta[im]  ->Fill(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it));
	  h_resVsPt_phi[im]  ->Fill(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it));
	  h_resVsPt_z0[im]   ->Fill(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it));

	  if (fabs(matchID_trk_eta->at(it)) < 0.8) {
	    h_resVsPt_pt_C[im]   ->Fill(matchID_trk_pt->at(it)  - matchID_simtrk_pt->at(it));
	    h_resVsPt_ptRel_C[im]->Fill((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it));
	    h_resVsPt_z0_C[im]   ->Fill(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it));
	    h_resVsPt_phi_C[im]  ->Fill(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it));
	  }
	  else if (fabs(matchID_trk_eta->at(it)) < 1.6 && fabs(matchID_trk_eta->at(it)) >= 0.8) {
	    h_resVsPt_pt_I[im]   ->Fill(matchID_trk_pt->at(it)  - matchID_simtrk_pt->at(it));
	    h_resVsPt_ptRel_I[im]->Fill((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it));
	    h_resVsPt_z0_I[im]   ->Fill(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it));
	    h_resVsPt_phi_I[im]  ->Fill(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it));
	  }
	  else if (fabs(matchID_trk_eta->at(it)) >= 1.6) {
	    h_resVsPt_pt_F[im]   ->Fill(matchID_trk_pt->at(it)  - matchID_simtrk_pt->at(it));
	    h_resVsPt_ptRel_F[im]->Fill((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it));
	    h_resVsPt_z0_F[im]   ->Fill(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it));
	    h_resVsPt_phi_F[im]  ->Fill(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it));
	  }

	}
      }

      // fill resolution vs. eta histograms
      for (int im=0; im<nETARANGE; im++) {
       if ( (fabs(matchID_simtrk_eta->at(it)) > (float)im*0.1) && (fabs(matchID_simtrk_eta->at(it)) < (float)im*0.1+0.1) ) {
	 h_resVsEta_pt[im]   ->Fill(matchID_trk_pt->at(it)  - matchID_simtrk_pt->at(it));
	 h_resVsEta_ptRel[im]->Fill((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it));
	 h_resVsEta_eta[im]  ->Fill(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it));
	 h_resVsEta_phi[im]  ->Fill(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it));
	 h_resVsEta_z0[im]   ->Fill(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it));

	 if (matchID_trk_pt->at(it)<5) h_resVsEta_eta_L[im]->Fill(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it));
	 else if (matchID_trk_pt->at(it)<15 && matchID_trk_pt->at(it)>=5) h_resVsEta_eta_M[im]->Fill(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it));
	 else h_resVsEta_eta_H[im]->Fill(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it));

       }
      }
      
      // fill 2D histograms
      h_2d_dz0_eta   ->Fill(matchID_simtrk_eta->at(it), fabs(matchID_trk_z0->at(it)  - matchID_simtrk_z0->at(it)));
      h_2d_deta_eta  ->Fill(matchID_simtrk_eta->at(it), fabs(matchID_trk_eta->at(it) - matchID_simtrk_eta->at(it)));
      h_2d_dphi_eta  ->Fill(matchID_simtrk_eta->at(it), fabs(matchID_trk_phi->at(it) - matchID_simtrk_phi->at(it)));
      h_2d_dptRel_eta->Fill(matchID_simtrk_eta->at(it), fabs((matchID_trk_pt->at(it) - matchID_simtrk_pt->at(it))/matchID_simtrk_pt->at(it)));
     
    } // end of matched track loop
  
  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------
  

  // ----------------------------------------------------------------------------------------------------------------
  // 2D plots  
  // ----------------------------------------------------------------------------------------------------------------

  TH1F* h2_resVsPt_pt   = new TH1F("resVsPt2_pt",   ";Sim track p_{T} [GeV]; p_{T} resolution", 20,0,100);
  TH1F* h2_resVsPt_pt_C = new TH1F("resVsPt2_pt_C", ";Sim track p_{T} [GeV]; p_{T} resolution", 20,0,100);
  TH1F* h2_resVsPt_pt_I = new TH1F("resVsPt2_pt_I", ";Sim track p_{T} [GeV]; p_{T} resolution", 20,0,100);
  TH1F* h2_resVsPt_pt_F = new TH1F("resVsPt2_pt_F", ";Sim track p_{T} [GeV]; p_{T} resolution", 20,0,100);

  TH1F* h2_resVsPt_ptRel = new TH1F("resVsPt2_ptRel",     ";Sim track p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_ptRel_C = new TH1F("resVsPt2_ptRel_C", ";Sim track p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_ptRel_I = new TH1F("resVsPt2_ptRel_I", ";Sim track p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);
  TH1F* h2_resVsPt_ptRel_F = new TH1F("resVsPt2_ptRel_F", ";Sim track p_{T} [GeV]; p_{T} resolution / p_{T}", 20,0,100);

  TH1F* h2_mresVsPt_pt   = new TH1F("mresVsPt2_pt",   ";Sim track p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);
  TH1F* h2_mresVsPt_pt_C = new TH1F("mresVsPt2_pt_C", ";Sim track p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);
  TH1F* h2_mresVsPt_pt_I = new TH1F("mresVsPt2_pt_I", ";Sim track p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);
  TH1F* h2_mresVsPt_pt_F = new TH1F("mresVsPt2_pt_F", ";Sim track p_{T} [GeV]; Mean(p_{T} residual) [GeV]", 20,0,100);

  TH1F* h2_resVsPt_z0   = new TH1F("resVsPt2_z0",   ";Sim track p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);
  TH1F* h2_resVsPt_z0_C = new TH1F("resVsPt2_z0_C", ";Sim track p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);
  TH1F* h2_resVsPt_z0_I = new TH1F("resVsPt2_z0_I", ";Sim track p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);
  TH1F* h2_resVsPt_z0_F = new TH1F("resVsPt2_z0_F", ";Sim track p_{T} [GeV]; z_{0} resolution [cm]", 20,0,100);

  TH1F* h2_resVsPt_phi   = new TH1F("resVsPt2_phi",   ";Sim track p_{T} [GeV]; #phi resolution [rad]", 20,0,100);
  TH1F* h2_resVsPt_phi_C = new TH1F("resVsPt2_phi_C", ";Sim track p_{T} [GeV]; #phi resolution [rad]", 20,0,100);
  TH1F* h2_resVsPt_phi_I = new TH1F("resVsPt2_phi_I", ";Sim track p_{T} [GeV]; #phi resolution [rad]", 20,0,100);
  TH1F* h2_resVsPt_phi_F = new TH1F("resVsPt2_phi_F", ";Sim track p_{T} [GeV]; #phi resolution [rad]", 20,0,100);

  TH1F* h2_resVsPt_eta   = new TH1F("resVsPt2_eta",  ";Sim track p_{T} [GeV]; #eta resolution", 20,0,100);

  for (int i=0; i<nRANGE; i++) {
    // set bin content and error
    h2_resVsPt_pt  ->SetBinContent(i+1, h_resVsPt_pt[i]  ->GetRMS());
    h2_resVsPt_pt  ->SetBinError(  i+1, h_resVsPt_pt[i]  ->GetRMSError());
    h2_resVsPt_pt_C->SetBinContent(i+1, h_resVsPt_pt_C[i]->GetRMS());
    h2_resVsPt_pt_C->SetBinError(  i+1, h_resVsPt_pt_C[i]->GetRMSError());
    h2_resVsPt_pt_I->SetBinContent(i+1, h_resVsPt_pt_I[i]->GetRMS());
    h2_resVsPt_pt_I->SetBinError(  i+1, h_resVsPt_pt_I[i]->GetRMSError());
    h2_resVsPt_pt_F->SetBinContent(i+1, h_resVsPt_pt_F[i]->GetRMS());
    h2_resVsPt_pt_F->SetBinError(  i+1, h_resVsPt_pt_F[i]->GetRMSError());

    h2_resVsPt_ptRel  ->SetBinContent(i+1, h_resVsPt_ptRel[i]  ->GetRMS());
    h2_resVsPt_ptRel  ->SetBinError(  i+1, h_resVsPt_ptRel[i]  ->GetRMSError());
    h2_resVsPt_ptRel_C->SetBinContent(i+1, h_resVsPt_ptRel_C[i]->GetRMS());
    h2_resVsPt_ptRel_C->SetBinError(  i+1, h_resVsPt_ptRel_C[i]->GetRMSError());
    h2_resVsPt_ptRel_I->SetBinContent(i+1, h_resVsPt_ptRel_I[i]->GetRMS());
    h2_resVsPt_ptRel_I->SetBinError(  i+1, h_resVsPt_ptRel_I[i]->GetRMSError());
    h2_resVsPt_ptRel_F->SetBinContent(i+1, h_resVsPt_ptRel_F[i]->GetRMS());
    h2_resVsPt_ptRel_F->SetBinError(  i+1, h_resVsPt_ptRel_F[i]->GetRMSError());

    h2_mresVsPt_pt  ->SetBinContent(i+1, h_resVsPt_pt[i]  ->GetMean());
    h2_mresVsPt_pt  ->SetBinError(  i+1, h_resVsPt_pt[i]  ->GetMeanError());
    h2_mresVsPt_pt_C->SetBinContent(i+1, h_resVsPt_pt_C[i]->GetMean());
    h2_mresVsPt_pt_C->SetBinError(  i+1, h_resVsPt_pt_C[i]->GetMeanError());
    h2_mresVsPt_pt_I->SetBinContent(i+1, h_resVsPt_pt_I[i]->GetMean());
    h2_mresVsPt_pt_I->SetBinError(  i+1, h_resVsPt_pt_I[i]->GetMeanError());
    h2_mresVsPt_pt_F->SetBinContent(i+1, h_resVsPt_pt_F[i]->GetMean());
    h2_mresVsPt_pt_F->SetBinError(  i+1, h_resVsPt_pt_F[i]->GetMeanError());

    h2_resVsPt_z0  ->SetBinContent(i+1, h_resVsPt_z0[i]  ->GetRMS());
    h2_resVsPt_z0  ->SetBinError(  i+1, h_resVsPt_z0[i]  ->GetRMSError());
    h2_resVsPt_z0_C->SetBinContent(i+1, h_resVsPt_z0_C[i]->GetRMS());
    h2_resVsPt_z0_C->SetBinError(  i+1, h_resVsPt_z0_C[i]->GetRMSError());
    h2_resVsPt_z0_I->SetBinContent(i+1, h_resVsPt_z0_I[i]->GetRMS());
    h2_resVsPt_z0_I->SetBinError(  i+1, h_resVsPt_z0_I[i]->GetRMSError());
    h2_resVsPt_z0_F->SetBinContent(i+1, h_resVsPt_z0_F[i]->GetRMS());
    h2_resVsPt_z0_F->SetBinError(  i+1, h_resVsPt_z0_F[i]->GetRMSError());

    h2_resVsPt_phi  ->SetBinContent(i+1, h_resVsPt_phi[i]  ->GetRMS());
    h2_resVsPt_phi  ->SetBinError(  i+1, h_resVsPt_phi[i]  ->GetRMSError());
    h2_resVsPt_phi_C->SetBinContent(i+1, h_resVsPt_phi_C[i]->GetRMS());
    h2_resVsPt_phi_C->SetBinError(  i+1, h_resVsPt_phi_C[i]->GetRMSError());
    h2_resVsPt_phi_I->SetBinContent(i+1, h_resVsPt_phi_I[i]->GetRMS());
    h2_resVsPt_phi_I->SetBinError(  i+1, h_resVsPt_phi_I[i]->GetRMSError());
    h2_resVsPt_phi_F->SetBinContent(i+1, h_resVsPt_phi_F[i]->GetRMS());
    h2_resVsPt_phi_F->SetBinError(  i+1, h_resVsPt_phi_F[i]->GetRMSError());

    h2_resVsPt_eta->SetBinContent(i+1, h_resVsPt_eta[i]->GetRMS());
    h2_resVsPt_eta->SetBinError(  i+1, h_resVsPt_eta[i]->GetRMSError());
  }


  // resolution vs. eta histograms
  TH1F* h2_resVsEta_eta   = new TH1F("resVsEta2_eta",   ";Sim track |#eta|; #eta resolution", 25,0,2.5);
  TH1F* h2_resVsEta_eta_L = new TH1F("resVsEta2_eta_L", ";Sim track |#eta|; #eta resolution", 25,0,2.5);
  TH1F* h2_resVsEta_eta_M = new TH1F("resVsEta2_eta_M", ";Sim track |#eta|; #eta resolution", 25,0,2.5);
  TH1F* h2_resVsEta_eta_H = new TH1F("resVsEta2_eta_H", ";Sim track |#eta|; #eta resolution", 25,0,2.5);

  TH1F* h2_mresVsEta_eta   = new TH1F("mresVsEta2_eta",   ";Sim track |#eta|; Mean(#eta residual)", 25,0,2.5);
  TH1F* h2_mresVsEta_eta_L = new TH1F("mresVsEta2_eta_L", ";Sim track |#eta|; Mean(#eta residual)", 25,0,2.5);
  TH1F* h2_mresVsEta_eta_M = new TH1F("mresVsEta2_eta_M", ";Sim track |#eta|; Mean(#eta residual)", 25,0,2.5);
  TH1F* h2_mresVsEta_eta_H = new TH1F("mresVsEta2_eta_H", ";Sim track |#eta|; Mean(#eta residual)", 25,0,2.5);

  TH1F* h2_resVsEta_z0  = new TH1F("resVsEta2_z0",  ";Sim track |#eta|; z_{0} resolution [cm]", 25,0,2.5);
  TH1F* h2_resVsEta_phi = new TH1F("resVsEta2_phi", ";Sim track |#eta|; #phi resolution [rad]", 25,0,2.5);

  TH1F* h2_resVsEta_pt    = new TH1F("resVsEta2_pt",    ";Sim track |#eta|; p_{T} resolution [GeV]",   25,0,2.5);
  TH1F* h2_resVsEta_ptRel = new TH1F("resVsEta2_ptRel", ";Sim track |#eta|; p_{T} resolution / p_{T}", 25,0,2.5);

  for (int i=0; i<nETARANGE; i++) {
    // set bin content and error
    h2_resVsEta_eta  ->SetBinContent(i+1, h_resVsEta_eta[i]  ->GetRMS());
    h2_resVsEta_eta  ->SetBinError(  i+1, h_resVsEta_eta[i]  ->GetRMSError());
    h2_resVsEta_eta_L->SetBinContent(i+1, h_resVsEta_eta_L[i]->GetRMS());
    h2_resVsEta_eta_L->SetBinError(  i+1, h_resVsEta_eta_L[i]->GetRMSError());
    h2_resVsEta_eta_M->SetBinContent(i+1, h_resVsEta_eta_M[i]->GetRMS());
    h2_resVsEta_eta_M->SetBinError(  i+1, h_resVsEta_eta_M[i]->GetRMSError());
    h2_resVsEta_eta_H->SetBinContent(i+1, h_resVsEta_eta_H[i]->GetRMS());
    h2_resVsEta_eta_H->SetBinError(  i+1, h_resVsEta_eta_H[i]->GetRMSError());

    h2_mresVsEta_eta  ->SetBinContent(i+1, h_resVsEta_eta[i]  ->GetMean());
    h2_mresVsEta_eta  ->SetBinError(  i+1, h_resVsEta_eta[i]  ->GetMeanError());
    h2_mresVsEta_eta_L->SetBinContent(i+1, h_resVsEta_eta_L[i]->GetMean());
    h2_mresVsEta_eta_L->SetBinError(  i+1, h_resVsEta_eta_L[i]->GetMeanError());
    h2_mresVsEta_eta_M->SetBinContent(i+1, h_resVsEta_eta_M[i]->GetMean());
    h2_mresVsEta_eta_M->SetBinError(  i+1, h_resVsEta_eta_M[i]->GetMeanError());
    h2_mresVsEta_eta_H->SetBinContent(i+1, h_resVsEta_eta_H[i]->GetMean());
    h2_mresVsEta_eta_H->SetBinError(  i+1, h_resVsEta_eta_H[i]->GetMeanError());

    h2_resVsEta_z0 ->SetBinContent(i+1, h_resVsEta_z0[i] ->GetRMS());
    h2_resVsEta_z0 ->SetBinError(  i+1, h_resVsEta_z0[i] ->GetRMSError());
    h2_resVsEta_phi->SetBinContent(i+1, h_resVsEta_phi[i]->GetRMS());
    h2_resVsEta_phi->SetBinError(  i+1, h_resVsEta_phi[i]->GetRMSError());

    h2_resVsEta_pt   ->SetBinContent(i+1, h_resVsEta_pt[i]   ->GetRMS());
    h2_resVsEta_pt   ->SetBinError(  i+1, h_resVsEta_pt[i]   ->GetRMSError());
    h2_resVsEta_ptRel->SetBinContent(i+1, h_resVsEta_ptRel[i]->GetRMS());
    h2_resVsEta_ptRel->SetBinError(  i+1, h_resVsEta_ptRel[i]->GetRMSError());
  }


  // set minimum to zero
  h2_resVsPt_pt  ->SetMinimum(0);
  h2_resVsPt_pt_C->SetMinimum(0);
  h2_resVsPt_pt_I->SetMinimum(0);
  h2_resVsPt_pt_F->SetMinimum(0);

  h2_resVsPt_ptRel  ->SetMinimum(0);
  h2_resVsPt_ptRel_C->SetMinimum(0);
  h2_resVsPt_ptRel_I->SetMinimum(0);
  h2_resVsPt_ptRel_F->SetMinimum(0);

  h2_resVsPt_z0  ->SetMinimum(0);
  h2_resVsPt_z0_C->SetMinimum(0);
  h2_resVsPt_z0_I->SetMinimum(0);
  h2_resVsPt_z0_F->SetMinimum(0);

  h2_resVsPt_phi  ->SetMinimum(0);
  h2_resVsPt_phi_C->SetMinimum(0);
  h2_resVsPt_phi_I->SetMinimum(0);
  h2_resVsPt_phi_F->SetMinimum(0);

  h2_resVsPt_eta->SetMinimum(0);

  h2_resVsEta_eta  ->SetMinimum(0);
  h2_resVsEta_eta_L->SetMinimum(0);
  h2_resVsEta_eta_M->SetMinimum(0);
  h2_resVsEta_eta_H->SetMinimum(0);

  h2_resVsEta_z0 ->SetMinimum(0);
  h2_resVsEta_phi->SetMinimum(0);
  h2_resVsEta_pt   ->SetMinimum(0);
  h2_resVsEta_ptRel->SetMinimum(0);


  // draw and save plots
  char ctxt[500];
  TCanvas c("c");

  h2_resVsPt_pt->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsPt_pt.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsPt_pt.png");
  c.SaveAs(type+"_canvas.pdf("); // keep the pdf file open

  if (doDetailedPlots) {
    h2_resVsPt_pt_C->Draw();
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_pt_C.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_pt_C.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_pt_I->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_pt_I.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_pt_I.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_pt_F->Draw();
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_pt_F.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_pt_F.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h2_resVsPt_ptRel->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel.png");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h2_resVsPt_ptRel_C->Draw();
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel_C.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel_C.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_ptRel_I->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel_I.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel_I.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_ptRel_F->Draw();
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel_F.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_ptRel_F.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h2_mresVsPt_pt->Draw();
  c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt.eps");
  c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt.png");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h2_mresVsPt_pt_C->Draw();
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_C.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_C.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_mresVsPt_pt_I->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_I.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_I.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_mresVsPt_pt_F->Draw();
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_F.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_F.png");
    c.SaveAs(type+"_canvas.pdf");
    
    // combined plot for the pt residual bias
    h2_mresVsPt_pt_C->GetYaxis()->SetRangeUser(-1.5,1.5);
    h2_mresVsPt_pt_C->Draw();
    h2_mresVsPt_pt_I->Draw("same");
    h2_mresVsPt_pt_I->SetMarkerColor(kRed);
    h2_mresVsPt_pt_I->SetLineColor(kRed);
    h2_mresVsPt_pt_F->Draw("same");
    h2_mresVsPt_pt_F->SetMarkerColor(kBlue);
    h2_mresVsPt_pt_F->SetLineColor(kBlue);
    
    TLegend* l = new TLegend(0.22,0.22,0.42,0.42);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextSize(0.04);
    l->AddEntry(h2_mresVsPt_pt_C," |eta| < 0.8","l");
    l->AddEntry(h2_mresVsPt_pt_I," 0.8 < |eta| < 1.6","l");
    l->AddEntry(h2_mresVsPt_pt_F," |eta| > 1.6","l");
    l->SetTextFont(42);
    l->Draw();	

    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsPt_pt_comb.png");
    c.SaveAs(type+"_canvas.pdf");
  }
  
  h2_resVsPt_z0->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsPt_z0.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsPt_z0.png");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h2_resVsPt_z0_C->Draw();
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.22,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_z0_C.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_z0_C.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_z0_I->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.22,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_z0_I.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_z0_I.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_z0_F->Draw();
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.22,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_z0_F.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_z0_F.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h2_resVsPt_phi->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsPt_phi.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsPt_phi.png");
  c.SaveAs(type+"_canvas.pdf");


  if (doDetailedPlots) {
    h2_resVsPt_phi_C->Draw();
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_phi_C.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_phi_C.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_phi_I->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_phi_I.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_phi_I.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsPt_phi_F->Draw();
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsPt_phi_F.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsPt_phi_F.png");
    c.SaveAs(type+"_canvas.pdf");
  }
    
  h2_resVsPt_eta->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsPt_eta.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsPt_eta.png");
  c.SaveAs(type+"_canvas.pdf");


  h2_resVsEta_eta->Draw();
  sprintf(ctxt,"Full p_{T} range");
  mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs("TrkPlots/"+type+"_resVsEta_eta.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsEta_eta.png");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h2_resVsEta_eta_L->Draw();
    sprintf(ctxt,"p_{T} < 5 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsEta_eta_L.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsEta_eta_L.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsEta_eta_M->Draw();
    sprintf(ctxt,"5 < p_{T} < 15 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsEta_eta_M.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsEta_eta_M.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_resVsEta_eta_H->Draw();
    sprintf(ctxt,"p_{T} > 15 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_resVsEta_eta_H.eps");
    c.SaveAs("TrkPlots/"+type+"_resVsEta_eta_H.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h2_mresVsEta_eta->Draw();
  sprintf(ctxt,"Full p_{T} range");
  mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta.eps");
  c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta.png");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h2_mresVsEta_eta_L->Draw();
    sprintf(ctxt,"p_{T} < 5 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta_L.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta_L.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_mresVsEta_eta_M->Draw();
    sprintf(ctxt,"5 < p_{T} < 15 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta_M.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta_M.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h2_mresVsEta_eta_H->Draw();
    sprintf(ctxt,"p_{T} > 15 GeV");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta_H.eps");
    c.SaveAs("TrkPlots/"+type+"_mresVsEta_eta_H.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h2_resVsEta_z0->Draw();
  sprintf(ctxt,"Full p_{T} range");
  mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs("TrkPlots/"+type+"_resVsEta_z0.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsEta_z0.png");
  c.SaveAs(type+"_canvas.pdf");

  h2_resVsEta_phi->Draw();
  sprintf(ctxt,"Full p_{T} range");
  mySmallText(0.22,0.82,1,ctxt);
  c.SaveAs("TrkPlots/"+type+"_resVsEta_phi.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsEta_phi.png");
  c.SaveAs(type+"_canvas.pdf");

  h2_resVsEta_pt->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsEta_pt.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsEta_pt.png");
  c.SaveAs(type+"_canvas.pdf");

  h2_resVsEta_ptRel->Draw();
  c.SaveAs("TrkPlots/"+type+"_resVsEta_ptRel.eps");
  c.SaveAs("TrkPlots/"+type+"_resVsEta_ptRel.png");
  c.SaveAs(type+"_canvas.pdf");

  
  // ----------------------------------------------------------------------------------------------------------------
  // track quality plots
  // ----------------------------------------------------------------------------------------------------------------


  if (doDetailedPlots) {

    // draw and save plots
    h_matchID_trk_nstub->Draw();
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_nstub_C->Draw();
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub_C.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub_C.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_nstub_I->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub_I.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub_I.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_nstub_F->Draw();
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub_F.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_nstub_F.eps");
    c.SaveAs(type+"_canvas.pdf");
  }
    
  h_matchID_trk_chi2->Draw();
  sprintf(ctxt,"|eta| < 2.5; Full p_{T} range");
  mySmallText(0.52,0.82,1,ctxt);
  c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2.png");
  c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2.eps");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h_matchID_trk_chi2_C_L->Draw();
    sprintf(ctxt,"|eta| < 0.8; p_{T} < 5 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_L.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_L.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_I_L->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6; p_{T} < 5 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_L.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_L.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_F_L->Draw();
    sprintf(ctxt,"|eta| > 1.6; p_{T} < 5 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_L.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_L.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_C_M->Draw();
    sprintf(ctxt,"|eta| < 0.8; 5 < p_{T} < 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_M.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_M.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_I_M->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6, 5 < p_{T} < 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_M.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_M.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_F_M->Draw();
    sprintf(ctxt,"|eta| > 1.6; 5 < p_{T} < 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_M.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_M.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_C_H->Draw();
    sprintf(ctxt,"|eta| < 0.8; p_{T} > 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_H.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_H.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_I_H->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6; p_{T} > 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_H.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_H.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_F_H->Draw();
    sprintf(ctxt,"|eta| > 1.6; p_{T} > 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_H.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_H.eps");
    c.SaveAs(type+"_canvas.pdf");
  }
    
  h_matchID_trk_chi2_dof->Draw();
  sprintf(ctxt,"|eta| < 2.5; Full p_{T} range");
  mySmallText(0.52,0.82,1,ctxt);
  c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof.png");
  c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof.eps");
  c.SaveAs(type+"_canvas.pdf");

  if (doDetailedPlots) {
    h_matchID_trk_chi2_dof_C_L->Draw();
    sprintf(ctxt,"|eta| < 0.8; p_{T} < 5 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_L.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_L.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_I_L->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6; p_{T} < 5 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_L.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_L.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_F_L->Draw();
    sprintf(ctxt,"|eta| > 1.6; p_{T} < 5 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_L.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_L.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_C_M->Draw();
    sprintf(ctxt,"|eta| < 0.8; 5 < p_{T} < 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_M.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_M.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_I_M->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6, 5 < p_{T} < 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_M.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_M.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_F_M->Draw();
    sprintf(ctxt,"|eta| > 1.6; 5 < p_{T} < 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_M.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_M.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_C_H->Draw();
    sprintf(ctxt,"|eta| < 0.8; p_{T} > 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_H.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_H.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_I_H->Draw();
    sprintf(ctxt,"0.8 < |eta| < 1.6; p_{T} > 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_H.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_H.eps");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_F_H->Draw();
    sprintf(ctxt,"|eta| > 1.6; p_{T} > 15 GeV");
    mySmallText(0.52,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_H.png");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_H.eps");
    c.SaveAs(type+"_canvas.pdf");
  }
   
 
  // ----------------------------------------------------------------------------------------------------------------
  // efficiency plots  
  // ----------------------------------------------------------------------------------------------------------------

  // rebin pt/phi plots
  h_simtrk_pt->Rebin(2);
  h_simtrk_phi->Rebin(2);
  h_matchID_simtrk_pt->Rebin(2);
  h_matchID_simtrk_phi->Rebin(2);


  // calculate the effeciency
  h_matchID_simtrk_pt->Sumw2();
  h_simtrk_pt->Sumw2();
  TH1F* h_eff_pt = (TH1F*) h_matchID_simtrk_pt->Clone();
  h_eff_pt->SetName("eff_pt");
  h_eff_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt->Divide(h_matchID_simtrk_pt, h_simtrk_pt, 1.0, 1.0, "B");

  h_matchID_simtrk_pt_L->Sumw2();
  h_simtrk_pt_L->Sumw2();
  TH1F* h_eff_pt_L = (TH1F*) h_matchID_simtrk_pt_L->Clone();
  h_eff_pt_L->SetName("eff_pt_L");
  h_eff_pt_L->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_L->Divide(h_matchID_simtrk_pt_L, h_simtrk_pt_L, 1.0, 1.0, "B");

  h_matchID_simtrk_eta->Sumw2();
  h_simtrk_eta->Sumw2();
  TH1F* h_eff_eta = (TH1F*) h_matchID_simtrk_eta->Clone();
  h_eff_eta->SetName("eff_eta");
  h_eff_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta->Divide(h_matchID_simtrk_eta, h_simtrk_eta, 1.0, 1.0, "B");
  
  h_matchID_simtrk_phi->Sumw2();
  h_simtrk_phi->Sumw2();
  TH1F* h_eff_phi = (TH1F*) h_matchID_simtrk_phi->Clone();
  h_eff_phi->SetName("eff_phi");
  h_eff_phi->GetYaxis()->SetTitle("Efficiency");
  h_eff_phi->Divide(h_matchID_simtrk_phi, h_simtrk_phi, 1.0, 1.0, "B");


  // set the axis range
  h_eff_pt  ->SetAxisRange(0,1.1,"Y");
  h_eff_pt_L->SetAxisRange(0,1.1,"Y");
  h_eff_eta ->SetAxisRange(0,1.1,"Y");
  h_eff_phi ->SetAxisRange(0,1.1,"Y");
  

  // draw and save plots
  h_eff_pt->Draw();
  c.SaveAs("TrkPlots/"+type+"_eff_pt.eps");
  c.SaveAs("TrkPlots/"+type+"_eff_pt.png");
  c.SaveAs(type+"_canvas.pdf");

  if (type.Contains("SingleMu")) {
    h_eff_pt->GetYaxis()->SetRangeUser(0.9,1.01); // zoomed-in plot
    c.SaveAs("TrkPlots/"+type+"_eff_pt_zoom.eps");
    c.SaveAs("TrkPlots/"+type+"_eff_pt_zoom.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h_eff_pt_L->Draw();
  c.SaveAs("TrkPlots/"+type+"_eff_pt_L.eps");
  c.SaveAs("TrkPlots/"+type+"_eff_pt_L.png");
  c.SaveAs(type+"_canvas.pdf");
  
  h_eff_eta->Draw();
  c.SaveAs("TrkPlots/"+type+"_eff_eta.eps");
  c.SaveAs("TrkPlots/"+type+"_eff_eta.png");
  c.SaveAs(type+"_canvas.pdf");
  
  if (type.Contains("SingleMu")) {
    h_eff_eta->GetYaxis()->SetRangeUser(0.8,1.01); // zoomed-in plot
    c.SaveAs("TrkPlots/"+type+"_eff_eta_zoom.eps");
    c.SaveAs("TrkPlots/"+type+"_eff_eta_zoom.png");
    c.SaveAs(type+"_canvas.pdf");
  }

  h_eff_phi->Draw();
  c.SaveAs("TrkPlots/"+type+"_eff_phi.eps");
  c.SaveAs("TrkPlots/"+type+"_eff_phi.png");
  c.SaveAs(type+"_canvas.pdf");
  
  if (type.Contains("SingleMu")) {
    h_eff_phi->GetYaxis()->SetRangeUser(0.8,1.01); // zoomed-in plot
    c.SaveAs("TrkPlots/"+type+"_eff_phi_zoom.eps");
    c.SaveAs("TrkPlots/"+type+"_eff_phi_zoom.png");
    c.SaveAs(type+"_canvas.pdf");
  }


  // ----------------------------------------------------------------------------------------------------------------
  // resolution plots
  // ----------------------------------------------------------------------------------------------------------------

  float rms = 0;

  if (doDetailedPlots) {

    // draw and save plots
    h_res_pt->Draw();
    rms = h_res_pt->GetRMS();
    sprintf(ctxt,"RMS = %.4f",rms);
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_pt.eps");
    c.SaveAs("TrkPlots/"+type+"_res_pt.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_ptRel->Draw();
    rms = h_res_ptRel->GetRMS();
    sprintf(ctxt,"RMS = %.4f",rms);	
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_ptRel.eps");
    c.SaveAs("TrkPlots/"+type+"_res_ptRel.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_eta->Draw();
    rms = h_res_eta->GetRMS();
    sprintf(ctxt,"RMS = %.3e",rms);	
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_res_eta.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_phi->Draw();
    rms = h_res_phi->GetRMS();
    sprintf(ctxt,"RMS = %.3e",rms);	
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_phi.eps");
    c.SaveAs("TrkPlots/"+type+"_res_phi.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_z0->Draw();
    rms = h_res_z0->GetRMS();
    sprintf(ctxt,"RMS = %.4f",rms);	
    mySmallText(0.22,0.82,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_z0.eps");
    c.SaveAs("TrkPlots/"+type+"_res_z0.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_z0_C->Draw();
    rms = h_res_z0_C->GetRMS();
    sprintf(ctxt,"RMS = %.4f;",rms);
    mySmallText(0.22,0.82,1,ctxt);
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_z0_C.eps");
    c.SaveAs("TrkPlots/"+type+"_res_z0_C.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_z0_I->Draw();
    rms = h_res_z0_I->GetRMS();
    sprintf(ctxt,"RMS = %.4f;",rms);
    mySmallText(0.22,0.82,1,ctxt);
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_z0_I.eps");
    c.SaveAs("TrkPlots/"+type+"_res_z0_I.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_res_z0_F->Draw();
    rms = h_res_z0_F->GetRMS();
    sprintf(ctxt,"RMS = %.4f;",rms);
    mySmallText(0.22,0.82,1,ctxt);
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.22,0.76,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_res_z0_F.eps");
    c.SaveAs("TrkPlots/"+type+"_res_z0_F.png");
    c.SaveAs(type+"_canvas.pdf");
  }
   
 
  // ----------------------------------------------------------------------------------------------------------------
  // 2D histogram plots
  // ----------------------------------------------------------------------------------------------------------------

  if (doDetailedPlots) {

    // draw and save plots
    h_2d_logchi2_eta->Draw("colz");
    c.SaveAs("TrkPlots/"+type+"_2d_logchi2_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_2d_logchi2_eta.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_2d_logchi2_dof_eta->Draw("colz");
    c.SaveAs("TrkPlots/"+type+"_2d_logchi2_dof_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_2d_logchi2_dof_eta.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_2d_dz0_eta->Draw("colz");
    c.SaveAs("TrkPlots/"+type+"_2d_dz0_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_2d_dz0_eta.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_2d_deta_eta->Draw("colz");
    c.SaveAs("TrkPlots/"+type+"_2d_deta_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_2d_deta_eta.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_2d_dphi_eta->Draw("colz");
    c.SaveAs("TrkPlots/"+type+"_2d_dphi_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_2d_dphi_eta.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_2d_dptRel_eta->Draw("colz");
    c.SaveAs("TrkPlots/"+type+"_2d_dptRel_eta.eps");
    c.SaveAs("TrkPlots/"+type+"_2d_dptRel_eta.png");
    c.SaveAs(type+"_canvas.pdf");
  }
    
  // ----------------------------------------------------------------------------------------------------------------
  // logarithmic plots
  // ----------------------------------------------------------------------------------------------------------------
  
  // set y-axis to a log scale
  c.SetLogy();
    
  if (doDetailedPlots) {
    
    // draw and save plots
    h_matchID_trk_chi2->Draw();
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_log.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_log.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_C_H->Draw();
    h_matchID_trk_chi2_C_M->Draw("same");
    h_matchID_trk_chi2_C_M->SetLineColor(kRed);
    h_matchID_trk_chi2_C_L->Draw("same");
    h_matchID_trk_chi2_C_L->SetLineColor(kBlue);

    TLegend* lh = new TLegend(0.65,0.7,0.85,0.9);
    lh->SetFillStyle(0);
    lh->SetBorderSize(0);
    lh->SetTextSize(0.04);
    lh->AddEntry(h_matchID_trk_chi2_C_H," p_{T} > 15 GeV","l");
    lh->AddEntry(h_matchID_trk_chi2_C_M," 5 < p_{T} < 15 GeV","l");
    lh->AddEntry(h_matchID_trk_chi2_C_L," p_{T} < 5 GeV","l");
    lh->SetTextFont(42);
    lh->Draw();	
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.66,0.6,1,ctxt);

    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_log_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_C_log_comb.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_I_H->Draw();
    h_matchID_trk_chi2_I_M->Draw("same");
    h_matchID_trk_chi2_I_M->SetLineColor(kRed);
    h_matchID_trk_chi2_I_L->Draw("same");
    h_matchID_trk_chi2_I_L->SetLineColor(kBlue);
    lh->Draw();	
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.66,0.6,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_log_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_I_log_comb.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_F_H->Draw();
    h_matchID_trk_chi2_F_M->Draw("same");
    h_matchID_trk_chi2_F_M->SetLineColor(kRed);
    h_matchID_trk_chi2_F_L->Draw("same");
    h_matchID_trk_chi2_F_L->SetLineColor(kBlue);
    lh->Draw();	
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.66,0.6,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_log_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_F_log_comb.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof->Draw();
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_log.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_log.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_C_H->Draw();
    h_matchID_trk_chi2_dof_C_M->Draw("same");
    h_matchID_trk_chi2_dof_C_M->SetLineColor(kRed);
    h_matchID_trk_chi2_dof_C_L->Draw("same");
    h_matchID_trk_chi2_dof_C_L->SetLineColor(kBlue);
    lh->Draw();	
    sprintf(ctxt,"|eta| < 0.8");
    mySmallText(0.66,0.6,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_log_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_C_log_comb.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_I_H->Draw();
    h_matchID_trk_chi2_dof_I_M->Draw("same");
    h_matchID_trk_chi2_dof_I_M->SetLineColor(kRed);
    h_matchID_trk_chi2_dof_I_L->Draw("same");
    h_matchID_trk_chi2_dof_I_L->SetLineColor(kBlue);
    lh->Draw();	
    sprintf(ctxt,"0.8 < |eta| < 1.6");
    mySmallText(0.66,0.6,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_log_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_I_log_comb.png");
    c.SaveAs(type+"_canvas.pdf");
    
    h_matchID_trk_chi2_dof_F_H->Draw();
    h_matchID_trk_chi2_dof_F_M->Draw("same");
    h_matchID_trk_chi2_dof_F_M->SetLineColor(kRed);
    h_matchID_trk_chi2_dof_F_L->Draw("same");
    h_matchID_trk_chi2_dof_F_L->SetLineColor(kBlue);
    lh->Draw();	
    sprintf(ctxt,"|eta| > 1.6");
    mySmallText(0.66,0.6,1,ctxt);
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_log_comb.eps");
    c.SaveAs("TrkPlots/"+type+"_matchID_trk_chi2_dof_F_log_comb.png");
    c.SaveAs(type+"_canvas.pdf)"); // close the pdf file
  }
  else {
    c.Clear();
    c.SaveAs(type+"_canvas.pdf)"); // close the pdf file
  }

}


void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


