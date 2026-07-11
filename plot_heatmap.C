#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include <cmath>
#include <iostream>

void plot_heatmap(int eventToShow = 0) {

    TFile *f = TFile::Open("full_history_pythia.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open full_history_pythia.root" << std::endl;
        return;
    }
    TTree *t = (TTree*)f->Get("ParticleTree");

    int EventID, ParticleID;
    bool IsFinal;
    double pT, eta, phi, Energy;
    
    t->SetBranchAddress("EventID", &EventID);
    t->SetBranchAddress("ParticleID", &ParticleID);
    t->SetBranchAddress("IsFinal", &IsFinal);
    t->SetBranchAddress("pT", &pT);
    t->SetBranchAddress("eta", &eta);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("Energy", &Energy);

    TH2F *h_Et = new TH2F("h_Et", ";#eta;#phi;E_T [GeV]", 40, -5.0, 5.0, 40, -M_PI, M_PI);
    TH2F *h_pT = new TH2F("h_pT", ";#eta;#phi;p_T [GeV]", 40, -5.0, 5.0, 40, -M_PI, M_PI);

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        
        if (EventID != eventToShow) continue;
        if (!IsFinal) continue; 
        
        int absPDG = std::abs(ParticleID);
        if (absPDG == 12 || absPDG == 14 || absPDG == 16) continue;

        double Et = Energy / std::cosh(eta);

        h_Et->Fill(eta, phi, Et);
        h_pT->Fill(eta, phi, pT);
    }

    gStyle->SetOptStat(0);         
    gStyle->SetOptTitle(0);       
    gStyle->SetPalette(kBird);     
    gStyle->SetNumberContours(256); 
    gStyle->SetCanvasColor(kWhite);

    //undeflow count and overflow count
    double et_under = 0, et_over = 0;
    double pt_under = 0, pt_over = 0;
    
    for (int j = 1; j <= h_Et->GetNbinsY(); j++) {
        et_under += h_Et->GetBinContent(0, j);                     
        et_over  += h_Et->GetBinContent(h_Et->GetNbinsX() + 1, j); 
        pt_under += h_pT->GetBinContent(0, j);                     
        pt_over  += h_pT->GetBinContent(h_pT->GetNbinsX() + 1, j); 
    }

   
    TCanvas *c1 = new TCanvas("c1", "Eta-Phi Heatmap", 1200, 600);
    c1->Divide(2, 1); 


    c1->cd(1);
    gPad->SetRightMargin(0.15); 
    gPad->SetTopMargin(0.15);  
    h_Et->Draw("COLZ");         

    TPaveText *pt_Et = new TPaveText(0.10, 0.86, 0.85, 0.98, "NDC");
    pt_Et->SetFillColor(kWhite);
    pt_Et->SetBorderSize(0);    
    pt_Et->SetTextAlign(12);    
    pt_Et->AddText(Form("Transverse Energy (E_T) - Event %d", eventToShow));
    pt_Et->AddText(Form("Underflow (#eta < -5): %.2f GeV  |  Overflow (#eta > 5): %.2f GeV", et_under, et_over));
    pt_Et->Draw();

   
    c1->cd(2);
    gPad->SetRightMargin(0.15); 
    gPad->SetTopMargin(0.15);
    h_pT->Draw("COLZ");

  
    TPaveText *pt_pT = new TPaveText(0.10, 0.86, 0.85, 0.98, "NDC");
    pt_pT->SetFillColor(kWhite);
    pt_pT->SetBorderSize(0);
    pt_pT->SetTextAlign(12);
    pt_pT->AddText(Form("Transverse Momentum (p_T) - Event %d", eventToShow));
    pt_pT->AddText(Form("Underflow (#eta < -5): %.2f GeV  |  Overflow (#eta > 5): %.2f GeV", pt_under, pt_over));
    pt_pT->Draw();

    
    c1->SaveAs(Form("Event_%d_Heatmap.png", eventToShow));
    c1->SaveAs(Form("Event_%d_Heatmap.pdf", eventToShow));
    
    std::cout << "Saved plots as Event_" << eventToShow << "_Heatmap.png and .pdf" << std::endl;
}