{

//TFile f("MinBias_Em_from_Run1_L1EG.root");

TChain* Events = new TChain("Events");

//Events -> Add("PROD/MinBias_L1TrackElectron.root");
//Events -> Add("PROD/MinBias_L1TrackElectron_FerdosSettings.root");
Events -> Add("PROD/MinBias_L1TrackElectron_FerdosSettingsLoose.root");


gStyle->SetOptStat(0);

int nbins=90;
float x1 = 4.5;
float x2 = 94.5;
	// Run-1 algos :
/*
TH1F* h1 = new TH1F("hIsoTtk",";ET (GeV); Events",nbins,x1,x2);
TH1F* h2 = new TH1F("hEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* h3 = new TH1F("hIsoEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* h4 = new TH1F("hIsoEGIsoTtk",";ET (GeV); Events",nbins,x1,x2);

TH1F* h1int = new TH1F("hIsoTtkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* h2int = new TH1F("hEGint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* h3int = new TH1F("hIsoEGint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* h4int = new TH1F("hIsoEGIsoTtkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);

Events->Draw("TMath::Max(TMath::Max(l1extraL1EmParticles_l1extraParticles_NonIsolated_Ele.obj.pt_[1],l1extraL1EmParticles_l1extraParticles_Isolated_Ele.obj.pt_[1]),  TMath::Max(l1extraL1EmParticles_l1extraParticles_NonIsolated_Ele.obj.pt_[0],l1extraL1EmParticles_l1extraParticles_Isolated_Ele.obj.pt_[0]) )>>hEG") ;

Events->Draw("TMath::Max(TMath::Max(l1extraL1TkElectronParticle_L1TkElectronsRun1EG_IsoTrk_Ele.obj.pt_[1],l1extraL1TkElectronParticle_L1TkElectronsRunIso1EG_IsoTrk_Ele.obj.pt_[1]), TMath::Max(l1extraL1TkElectronParticle_L1TkElectronsRun1EG_IsoTrk_Ele.obj.pt_[0],l1extraL1TkElectronParticle_L1TkElectronsRunIso1EG_IsoTrk_Ele.obj.pt_[0]))>>hIsoTtk") ;

Events->Draw("TMath::Max(l1extraL1EmParticles_l1extraParticles_Isolated_Ele.obj.pt_[1],l1extraL1EmParticles_l1extraParticles_Isolated_Ele.obj.pt_[0])>>hIsoEG") ;

Events->Draw("TMath::Max(l1extraL1TkElectronParticle_L1TkElectronsRunIso1EG_IsoTrk_Ele.obj.pt_[1],l1extraL1TkElectronParticle_L1TkElectronsRunIso1EG_IsoTrk_Ele.obj.pt_[0])>>hIsoEGIsoTtk") ;
*/

	// old Stage-2 algos :

TH1F* g1 = new TH1F("gEleTk",";ET (GeV); Events",nbins,x1,x2);
TH1F* g2 = new TH1F("gEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* g3 = new TH1F("gIsoEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* g4 = new TH1F("gIsoEGEleTk",";ET (GeV); Events",nbins,x1,x2);

TH1F* g5 = new TH1F("gEleTkIsoTk",";ET (GeV); Events",nbins,x1,x2);
TH1F* g6 = new TH1F("gIsoEGEleTkIsoTk",";ET (GeV); Events",nbins,x1,x2);


TH1F* g1int = new TH1F("gEleTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* g2int = new TH1F("gEGint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* g3int = new TH1F("gIsoEGint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* g4int = new TH1F("gIsoEGEleTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* g5int = new TH1F("gEleTkIsoTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* g6int = new TH1F("gIsoEGEleTkIsoTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);

/*
Events->Draw("TMath::Max(l1extraL1EmParticles_SLHCL1ExtraParticles_EGamma_Ele.obj.pt_[1],l1extraL1EmParticles_SLHCL1ExtraParticles_EGamma_Ele.obj.pt_[0])>>gEG") ;
Events->Draw("TMath::Max(l1extraL1TkElectronParticles_L1TkElectronsStage2EG__Ele.obj.pt_[1],l1extraL1TkElectronParticles_L1TkElectronsStage2EG__Ele.obj.pt_[0])>>gEleTk") ;
Events->Draw("TMath::Max(l1extraL1EmParticles_SLHCL1ExtraParticles_IsoEGamma_Ele.obj.pt_[1],l1extraL1EmParticles_SLHCL1ExtraParticles_IsoEGamma_Ele.obj.pt_[0])>>gIsoEG") ;
Events->Draw("TMath::Max(l1extraL1TkElectronParticles_L1TkElectronsStage2IsoEG__Ele.obj.pt_[1],l1extraL1TkElectronParticles_L1TkElectronsStage2IsoEG__Ele.obj.pt_[0])>>gIsoEGEleTk") ;

Events->Draw("TMath::Max(l1extraL1TkElectronParticles_L1TkElectronsStage2EGIsoTrk__Ele.obj.pt_[1],l1extraL1TkElectronParticles_L1TkElectronsStage2EGIsoTrk__Ele.obj.pt_[0])>>gEleTkIsoTk");
Events->Draw("TMath::Max(l1extraL1TkElectronParticles_L1TkElectronsStage2IsoEGIsoTrk__Ele.obj.pt_[1],l1extraL1TkElectronParticles_L1TkElectronsStage2IsoEGIsoTrk__Ele.obj.pt_[0])>>gIsoEGEleTkIsoTk");
*/

Events->Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticles_EGamma_Ele.obj.pt_)>>gEG");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkElectrons_EG_Ele.obj.pt_)>>gEleTk");
Events->Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticles_IsoEGamma_Ele.obj.pt_)>>gIsoEG");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkElectronsIsoEG_EG_Ele.obj.pt_)>>gIsoEGEleTk");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkIsoElectrons_EG_Ele.obj.pt_)>>gEleTkIsoTk");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkIsoElectronsIsoEG_EG_Ele.obj.pt_)>>gIsoEGEleTkIsoTk");


	// new clustering :

TH1F* j1 = new TH1F("jEleTk",";ET (GeV); Events",nbins,x1,x2);
TH1F* j2 = new TH1F("jEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* j3 = new TH1F("jIsoEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* j4 = new TH1F("jIsoEGEleTk",";ET (GeV); Events",nbins,x1,x2);
TH1F* j5 = new TH1F("jEleTkIsoTk",";ET (GeV); Events",nbins,x1,x2);
TH1F* j6 = new TH1F("jIsoEGEleTkIsoTk",";ET (GeV); Events",nbins,x1,x2);

TH1F* j1int = new TH1F("jEleTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* j2int = new TH1F("jEGint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* j3int = new TH1F("jIsoEGint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* j4int = new TH1F("jIsoEGEleTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* j5int = new TH1F("jEleTkIsoTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* j6int = new TH1F("jIsoEGEleTkIsoTkint",";ET threshold (GeV); Rate (kHz)",nbins,x1,x2);


Events->Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticlesNewClustering_EGamma_Ele.obj.pt_)>>jEG");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkElectronsNewclus_EG_Ele.obj.pt_)>>jEleTk");
Events->Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticlesNewClustering_IsoEGamma_Ele.obj.pt_)>>jIsoEG");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkElectronsNewclusIsoEG_EG_Ele.obj.pt_)>>jIsoEGEleTk");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkIsoElectronsNewclus_EG_Ele.obj.pt_)>>jEleTkIsoTk");
Events->Draw("Max$(l1extraL1TkElectronParticles_L1TkIsoElectronsNewclusIsoEG_EG_Ele.obj.pt_)>>jIsoEGEleTkIsoTk");


float nevts = 149500. ;
//float nevts = 20000;

/*
for (int i=0; i <= nbins+1; i++) {
  float v1 = h1->Integral(i,nbins+1); 
  float v2 = h2 -> Integral(i,nbins+1);
  float v3 = h3 -> Integral(i,nbins+1);
  float v4 = h4 -> Integral(i,nbins+1);

  h1int -> SetBinContent(i, v1);
  h2int -> SetBinContent(i, v2);
  h3int -> SetBinContent(i, v3);
  h4int -> SetBinContent(i, v4);
}

h1int -> Scale( 30000./nevts);
h2int -> Scale( 30000./nevts);
h3int -> Scale( 30000./nevts);
h4int -> Scale( 30000./nevts);
*/


for (int i=0; i <= nbins+1; i++) {
  float v1 = g1->Integral(i,nbins+1);
  float v2 = g2 -> Integral(i,nbins+1);
  float v3 = g3 -> Integral(i,nbins+1);
  float v4 = g4 -> Integral(i,nbins+1);
  float v5 = g5 -> Integral(i,nbins+1);
  float v6 = g6 -> Integral(i,nbins+1);

  g1int -> SetBinContent(i, v1);
  g2int -> SetBinContent(i, v2);
  g3int -> SetBinContent(i, v3);
  g4int -> SetBinContent(i, v4);
  g5int -> SetBinContent(i, v5);
  g6int -> SetBinContent(i, v6);

}

TH1F* ratio = new TH1F("ratio",";ET threshold (GeV); Rate reduction factor",nbins,x1,x2);
//ratio -> Divide(g1int, g2int, 1.0, 1.0, "B");
for (int i=0; i <= nbins; i++) {
  float den = g2int -> GetBinContent(i) ;
  float val =  ( g1int -> GetBinContent(i) ) / den ;
  float er = sqrt( val * (1.-val) / den );
  if (val > 0 && val < 1) {
     float rat = 1./ val ;
     float err = er / (val * val);
     ratio -> SetBinContent(i, rat);
     ratio -> SetBinError(i, err);
  }
  else {
     ratio -> SetBinContent(i, 0.);
  }
}


g1int -> Scale( 30000./nevts);
g2int -> Scale( 30000./nevts);
g3int -> Scale( 30000./nevts);
g4int -> Scale( 30000./nevts);
g5int -> Scale( 30000./nevts);
g6int -> Scale( 30000./nevts);


for (int i=0; i <= nbins+1; i++) {
  float v1 = j1->Integral(i,nbins+1);
  float v2 = j2 -> Integral(i,nbins+1);
  float v3 = j3 -> Integral(i,nbins+1);
  float v4 = j4 -> Integral(i,nbins+1);
  float v5 = j5 -> Integral(i,nbins+1);
  float v6 = j6 -> Integral(i,nbins+1);

  j1int -> SetBinContent(i, v1);
  j2int -> SetBinContent(i, v2);
  j3int -> SetBinContent(i, v3);
  j4int -> SetBinContent(i, v4);
  j5int -> SetBinContent(i, v5);
  j6int -> SetBinContent(i, v6);
}

j1int -> Scale( 30000./nevts);
j2int -> Scale( 30000./nevts);
j3int -> Scale( 30000./nevts);
j4int -> Scale( 30000./nevts);
j5int -> Scale( 30000./nevts); 
j6int -> Scale( 30000./nevts); 


float rmin = 10.;

/*
TCanvas* c1 = new TCanvas("c1","c1");

h1int -> SetMinimum(rmin);
h1int -> SetMaximum(40000.);
h1int -> Draw();
h2int -> SetLineColor(2);
h2int -> Draw("same");
h3int -> SetLineColor(4);
h3int -> Draw("same");
h4int -> SetLineColor(7);
h4int -> Draw("same");

gPad -> SetLogy(1);
gPad -> SetGridx(1);
gPad -> SetGridy(1);

 TLegend *leg = new TLegend(0.36,0.58,0.98,0.92,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
 leg -> SetHeader("Run 1 algorithms");
 leg -> AddEntry(g2int,"L1EG ","l");
 leg -> AddEntry(g3int,"L1IsoEG ","l");
 leg -> AddEntry(g1int,"L1TkElectron","l");
 leg -> AddEntry(g4int,"L1IsoEG_IsoTrk","l");
 leg -> Draw("same");
*/


TCanvas* c2 = new TCanvas("c2","c2");

g1int -> SetMinimum(rmin);
g1int -> SetMaximum(40000.);
g1int -> Draw();
g2int -> SetLineColor(2);
g2int -> Draw("same");
g3int -> SetLineColor(4);
g3int -> Draw("same");
g4int -> SetLineColor(7);
g4int -> Draw("same");
g5int -> SetLineColor(6);
g5int -> Draw("same");
g6int -> SetLineColor(3);
g6int -> Draw("same");

gPad -> SetLogy(1);
gPad -> SetGridx(1);
gPad -> SetGridy(1);

 TLegend *leg = new TLegend(0.36,0.58,0.98,0.92,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
 leg -> SetHeader("old stage-2 algorithms");
 leg -> AddEntry(g2int,"L1EG ","l");
 leg -> AddEntry(g3int,"L1IsoEG ","l");
 leg -> AddEntry(g1int,"L1TkElectron","l");
 leg -> AddEntry(g4int,"L1TkElectron_IsoEG","l");
 leg -> AddEntry(g5int,"L1TkElectron_IsoTk","l");
 leg -> AddEntry(g6int,"L1TkElectron_IsoEG_IsoTk","l");

 leg -> Draw("same");

        TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
        pt->SetBorderSize(0);
        pt->SetTextAlign(12);
        pt->SetFillStyle(0);
        pt->SetTextFont(42);
        pt->SetTextSize(0.03);
        //TText* text = pt->AddText(0.05,0.5,"Atanu's isolation settings");
        TText* text = pt->AddText(0.05,0.5,"Ferdos's isolation settings");
        pt->Draw();



TCanvas* c2bis = new TCanvas("c2bis","c2bis");
ratio -> Draw("pe");



TCanvas* c3 = new TCanvas("c3","c3");

j1int -> SetMinimum(rmin);
j1int -> SetMaximum(40000.);
j1int -> Draw();
j2int -> SetLineColor(2);
j2int -> Draw("same");
j3int -> SetLineColor(4);
j3int -> Draw("same");
j4int -> SetLineColor(7);
j4int -> Draw("same");
j5int -> SetLineColor(6);
j5int -> Draw("same");
j6int -> SetLineColor(3);
j6int -> Draw("same");

gPad -> SetLogy(1);
gPad -> SetGridx(1);
gPad -> SetGridy(1);

TLegend* leg3 = (TLegend*)leg->Clone();
leg3->SetHeader("new clustering algorithms");
leg3->Draw("same");
        pt->Draw();


}

