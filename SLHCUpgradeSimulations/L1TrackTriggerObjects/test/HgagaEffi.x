{

//TFile f("MinBias_Em_from_Run1_L1EG.root");

TChain* Events = new TChain("Events");
//Events -> Add("MinBias_Em_TrkMET_1.root");
//Events -> Add("MinBias_Em_TrkMET_2.root");
Events -> Add("Hgaga_Em_TrkMET_gen.root");

gStyle->SetOptStat(0);

int nbins=90;
float x1 = 4.5;
float x2 = 94.5;
	// Run-1 algos :

TH1F* h1 = new TH1F("hIsoTtk",";ET (GeV); Events",nbins,x1,x2);
TH1F* h2 = new TH1F("hEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* h3 = new TH1F("hIsoEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* h4 = new TH1F("hIsoEGIsoTtk",";ET (GeV); Events",nbins,x1,x2);

TH1F* h1int = new TH1F("hIsoTtkint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* h2int = new TH1F("hEGint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* h3int = new TH1F("hIsoEGint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* h4int = new TH1F("hIsoEGIsoTtkint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);

Events -> Draw("TMath::Max(Max$(l1extraL1EmParticles_l1extraParticles_NonIsolated_ALL.obj.pt_),Max$(l1extraL1EmParticles_l1extraParticles_Isolated_ALL.obj.pt_))>>hEG");

Events->Draw("TMath::Max(Max$(l1extraL1TkEmParticles_L1TkPhotonsRun1EG_IsoTrk_ALL.obj.pt_),Max$(l1extraL1TkEmParticles_L1TkPhotonsRunIso1EG_IsoTrk_ALL.obj.pt_))>>hIsoTtk");
Events -> Draw("Max$(l1extraL1EmParticles_l1extraParticles_Isolated_ALL.obj.pt_)>>hIsoEG");
Events -> Draw("Max$(l1extraL1TkEmParticles_L1TkPhotonsRunIso1EG_IsoTrk_ALL.obj.pt_)>>hIsoEGIsoTtk");


	// old Stage-2 algos :

TH1F* g1 = new TH1F("gIsoTtk",";ET (GeV); Events",nbins,x1,x2);
TH1F* g2 = new TH1F("gEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* g3 = new TH1F("gIsoEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* g4 = new TH1F("gIsoEGIsoTtk",";ET (GeV); Events",nbins,x1,x2);

TH1F* g1int = new TH1F("gIsoTtkint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* g2int = new TH1F("gEGint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* g3int = new TH1F("gIsoEGint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* g4int = new TH1F("gIsoEGIsoTtkint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);

Events -> Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticles_EGamma_ALL.obj.pt_)>>gEG");
Events -> Draw("Max$(l1extraL1TkEmParticles_L1TkPhotonsStage2EG_IsoTrk_ALL.obj.pt_)>>gIsoTtk");
Events -> Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticles_IsoEGamma_ALL.obj.pt_)>>gIsoEG");
Events -> Draw("Max$(l1extraL1TkEmParticles_L1TkPhotonsStage2IsoEG_IsoTrk_ALL.obj.pt_)>>gIsoEGIsoTtk");

	// new clustering :

TH1F* j1 = new TH1F("jIsoTtk",";ET (GeV); Events",nbins,x1,x2);
TH1F* j2 = new TH1F("jEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* j3 = new TH1F("jIsoEG",";ET (GeV); Events",nbins,x1,x2);
TH1F* j4 = new TH1F("jIsoEGIsoTtk",";ET (GeV); Events",nbins,x1,x2);

TH1F* j1int = new TH1F("jIsoTtkint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* j2int = new TH1F("jEGint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* j3int = new TH1F("jIsoEGint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);
TH1F* j4int = new TH1F("jIsoEGIsoTtkint",";ET threshold (GeV); SingleEG Efficiency",nbins,x1,x2);

Events -> Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticlesNewClustering_EGamma_ALL.obj.pt_)>>jEG");
Events -> Draw("Max$(l1extraL1TkEmParticles_L1TkPhotonsNewEG_IsoTrk_ALL.obj.pt_)>>jIsoTtk");
Events -> Draw("Max$(l1extraL1EmParticles_SLHCL1ExtraParticlesNewClustering_IsoEGamma_ALL.obj.pt_)>>jIsoEG");
Events -> Draw("Max$(l1extraL1TkEmParticles_L1TkPhotonsNewIsoEG_IsoTrk_ALL.obj.pt_)>>jIsoEGIsoTtk");



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

float nevts = 24099. ;

h1int -> Scale( 1./nevts);
h2int -> Scale( 1./nevts);
h3int -> Scale( 1./nevts);
h4int -> Scale( 1./nevts);


for (int i=0; i <= nbins+1; i++) {
  float v1 = g1->Integral(i,nbins+1);
  float v2 = g2 -> Integral(i,nbins+1);
  float v3 = g3 -> Integral(i,nbins+1);
  float v4 = g4 -> Integral(i,nbins+1);

  g1int -> SetBinContent(i, v1);
  g2int -> SetBinContent(i, v2);
  g3int -> SetBinContent(i, v3);
  g4int -> SetBinContent(i, v4);
}

g1int -> Scale( 1./nevts);
g2int -> Scale( 1./nevts);
g3int -> Scale( 1./nevts);
g4int -> Scale( 1./nevts);


for (int i=0; i <= nbins+1; i++) {
  float v1 = j1->Integral(i,nbins+1);
  float v2 = j2 -> Integral(i,nbins+1);
  float v3 = j3 -> Integral(i,nbins+1);
  float v4 = j4 -> Integral(i,nbins+1);

  j1int -> SetBinContent(i, v1);
  j2int -> SetBinContent(i, v2);
  j3int -> SetBinContent(i, v3);
  j4int -> SetBinContent(i, v4);
}

j1int -> Scale( 1./nevts);
j2int -> Scale( 1./nevts);
j3int -> Scale( 1./nevts);
j4int -> Scale( 1./nevts);


TCanvas* c1 = new TCanvas("c1","c1");

float rmin = 0.;
h1int -> SetMinimum(rmin);
h1int -> SetMaximum(1.1);
h1int -> Draw();
h2int -> SetLineColor(2);
h2int -> Draw("same");
h3int -> SetLineColor(4);
h3int -> Draw("same");
h4int -> SetLineColor(7);
h4int -> Draw("same");

//gPad -> SetLogy(1);
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
 leg -> AddEntry(h2int,"L1EG ","l");
 leg -> AddEntry(h3int,"L1IsoEG ","l");
 leg -> AddEntry(h1int,"L1EG_IsoTrk","l");
 leg -> AddEntry(h4int,"L1IsoEG_IsoTrk","l");
 leg -> Draw("same");

        TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
        pt->SetBorderSize(0);
        pt->SetTextAlign(12);
        pt->SetFillStyle(0);
        pt->SetTextFont(42);
        pt->SetTextSize(0.03);
        TText *text;
        text = pt->AddText(0.,0.2,"SingleEG Efficiency on Hgaga events");
        pt->Draw();



TCanvas* c2 = new TCanvas("c2","c2");

g1int -> SetMinimum(rmin);
g1int -> SetMaximum(1.1);
g1int -> Draw();
g2int -> SetLineColor(2);
g2int -> Draw("same");
g3int -> SetLineColor(4);
g3int -> Draw("same");
g4int -> SetLineColor(7);
g4int -> Draw("same");

//gPad -> SetLogy(1);
gPad -> SetGridx(1);
gPad -> SetGridy(1);
TLegend* leg2 = (TLegend*)leg->Clone();
leg2->SetHeader("old stage-2 algorithms");
leg2->Draw("same");
pt -> Draw();

TCanvas* c3 = new TCanvas("c3","c3");

j1int -> SetMinimum(rmin);
j1int -> SetMaximum(1.1);
j1int -> Draw();
j2int -> SetLineColor(2);
j2int -> Draw("same");
j3int -> SetLineColor(4);
j3int -> Draw("same");
j4int -> SetLineColor(7);
j4int -> Draw("same");

//gPad -> SetLogy(1);
gPad -> SetGridx(1);
gPad -> SetGridy(1);

TLegend* leg3 = (TLegend*)leg->Clone();
leg3->SetHeader("new clustering algorithms");
leg3->Draw("same");
pt -> Draw();

}

