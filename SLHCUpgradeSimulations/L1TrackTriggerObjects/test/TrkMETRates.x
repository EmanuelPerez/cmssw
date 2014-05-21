{


TChain* Events = new TChain("Events");
Events -> Add("MinBias_Em_TrkMET_1.root");
Events -> Add("MinBias_Em_TrkMET_2.root");

gStyle->SetOptStat(0);

int nbins=40;
float x1 = 0.;
float x2 = 100.;

TH1F* h1 = new TH1F("hTrkMET",";Missing ET (GeV); Events",nbins,x1,x2);
TH1F* h1int = new TH1F("hTrkMETint",";Missing ET threshold (GeV); Rate (kHz)",nbins,x1,x2);
TH1F* h2 = new TH1F("hMET",";Missing ET (GeV); Events",nbins,x1,x2);
TH1F* h2int = new TH1F("hMETint",";Missing ET threshold (GeV); Rate (kHz)",nbins,x1,x2);

Events -> Draw("l1extraL1TkEtMissParticles_L1TkEtMiss_MET_ALL.obj.pt_[0]>>hTrkMET");
Events -> Draw("l1extraL1EtMissParticles_l1extraParticles_MET_ALL.obj.pt_[0]>>hMET");


for (int i=0; i <= nbins+1; i++) {
  float v1 = h1->Integral(i,nbins+1); 
  h1int -> SetBinContent(i, v1);
  float v2 = h2->Integral(i,nbins+1); 
  h2int -> SetBinContent(i, v2);
}

float nevts = 149500. ;

h1int -> Scale( 30000./nevts);
h2int -> Scale( 30000./nevts);

TCanvas* c1 = new TCanvas("c1","c1");

float rmin = 10.;
h1int -> SetMinimum(rmin);
h1int -> SetMaximum(40000.);
h1int -> Draw();
h2int -> SetLineColor(2);
h2int->Draw("same");

gPad -> SetLogy(1);
gPad -> SetGridx(1);
gPad -> SetGridy(1);

 TLegend *leg = new TLegend(0.26,0.18,0.88,0.32,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
 leg -> AddEntry(h1int,"L1TrkMET rate","l");
 leg -> AddEntry(h2int,"L1MET rate","l");
 leg -> Draw("same");

TText t;
t.DrawTextNDC(0.2,0.35,"Scales are of course different !");



}

