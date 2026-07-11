#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>   // For writing to text files
#include <iomanip>   // For formatting the text file output

// FastJet Headers
#include "fastjet/ClusterSequence.hh"

// ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TEllipse.h"

using namespace std;
using namespace fastjet;

int main() {
    int eventToShow = 0;

    // 1. Open the ROOT file
    TFile *f = TFile::Open("full_history_pythia.root");
    if (!f || f->IsZombie()) {
        cerr << "Error: Could not open full_history_pythia.root" << endl;
        return 1;
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
    
    vector<PseudoJet> input_particles;

    // 2. Loop over the event
    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        
        if (EventID != eventToShow) continue;
        if (!IsFinal) continue; 
        
        // Skipping neutrinos
        int absPDG = abs(ParticleID);
        if (absPDG == 12 || absPDG == 14 || absPDG == 16) continue; // Skip neutrinos

        double Et = Energy / cosh(eta);
        h_Et->Fill(eta, phi, Et);
        h_pT->Fill(eta, phi, pT);

        double px = pT * cos(phi);
        double py = pT * sin(phi);
        double pz = pT * sinh(eta);
        
        PseudoJet particle(px, py, pz, Energy);
        particle.set_user_index(ParticleID); // Save the PDG ID inside the FastJet object
        input_particles.push_back(particle);
    }

    // 3. Run FastJet Clustering
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(input_particles, jet_def);
    
    double ptmin = 20;
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

    cout << "Found " << jets.size() << " jets with pT > " << ptmin << " GeV." << endl;

    // ========================================================================
    // 4. TEXT FILE OUTPUT: Write Jet and Constituent Data
    // ========================================================================
    string filename = "Event_" + to_string(eventToShow) + "_Jets_Info.txt";
    ofstream out(filename);
    
    out << "=========================================================\n";
    out << " Jet and Constituent Information for Event " << eventToShow << "\n";
    out << "=========================================================\n\n";

    for (unsigned i = 0; i < jets.size(); i++) {
        out << "---------------------------------------------------------\n";
        out << "JET [" << i << "]: "
            << "pT = " << fixed << setprecision(3) << jets[i].pt() << " GeV, "
            << "eta = " << jets[i].rap() << ", "
            << "phi = " << jets[i].phi_std() << "\n";
        out << "---------------------------------------------------------\n";

        vector<PseudoJet> constituents = jets[i].constituents();
        
        // Variables to manually calculate the sum of the constituents
        double sum_px = 0, sum_py = 0, sum_pz = 0, sum_E = 0;

        out << setw(10) << "PDG ID" << setw(15) << "pT [GeV]" 
            << setw(15) << "eta" << setw(15) << "phi" << "\n";
        
        for (unsigned j = 0; j < constituents.size(); j++) {
            int pdg_id = constituents[j].user_index();
            out << setw(10) << pdg_id
                << setw(15) << constituents[j].pt()
                << setw(15) << constituents[j].rap()
                << setw(15) << constituents[j].phi_std() << "\n";
                
            // Add up the 4-momentum components
            sum_px += constituents[j].px();
            sum_py += constituents[j].py();
            sum_pz += constituents[j].pz();
            sum_E  += constituents[j].E();
        }

        // Recreate the jet manually from the summed 4-momentum
        PseudoJet sum_jet(sum_px, sum_py, sum_pz, sum_E);
        
        out << ".........................................................\n";
        out << "SUM     : " 
            << "pT = " << sum_jet.pt() << " GeV, "
            << "eta = " << sum_jet.rap() << ", "
            << "phi = " << sum_jet.phi_std() << "\n";
        out << "---------------------------------------------------------\n\n";
    }
    
    out.close();
    cout << "Kinematics text file saved as " << filename << endl;

    // ========================================================================
    // 5. DRAWING THE HEATMAPS (Unchanged)
    // ========================================================================
    double et_under = 0, et_over = 0;
    double pt_under = 0, pt_over = 0;
    
    for (int j = 1; j <= h_Et->GetNbinsY(); j++) {
        et_under += h_Et->GetBinContent(0, j);                     
        et_over  += h_Et->GetBinContent(h_Et->GetNbinsX() + 1, j); 
        pt_under += h_pT->GetBinContent(0, j);                     
        pt_over  += h_pT->GetBinContent(h_pT->GetNbinsX() + 1, j); 
    }

    gStyle->SetOptStat(0);         
    gStyle->SetOptTitle(0);         
    gStyle->SetPalette(kBird);     
    gStyle->SetNumberContours(256); 
    gStyle->SetCanvasColor(kWhite);

    TCanvas *c1 = new TCanvas("c1", "Jet Heatmap", 1200, 600);
    c1->Divide(2, 1); 

    // --- Pad 1: Transverse Energy ---
    c1->cd(1);
    gPad->SetRightMargin(0.15); 
    gPad->SetTopMargin(0.15);   
    h_Et->Draw("COLZ");         

    TPaveText *pt_Et = new TPaveText(0.10, 0.86, 0.85, 0.98, "NDC");
    pt_Et->SetFillColor(kWhite);
    pt_Et->SetBorderSize(0);    
    pt_Et->SetTextAlign(12);    
    pt_Et->AddText(Form("Transverse Energy (E_T) - Event %d with R= %.1f and ptmin= %.1f GeV Jets", eventToShow, R, ptmin));
    pt_Et->AddText(Form("Underflow (#eta < -5): %.2f GeV  |  Overflow (#eta > 5): %.2f GeV", et_under, et_over));
    pt_Et->Draw();

    for (unsigned i = 0; i < jets.size(); i++) {
        TEllipse *circle = new TEllipse(jets[i].rap(), jets[i].phi_std(), R, R);
        circle->SetFillStyle(0); 
        circle->SetLineColor(kRed);
        circle->SetLineWidth(3);
        circle->Draw("same");
    }

    // --- Pad 2: Transverse Momentum ---
    c1->cd(2);
    gPad->SetRightMargin(0.15); 
    gPad->SetTopMargin(0.15);
    h_pT->Draw("COLZ");

    TPaveText *pt_pT = new TPaveText(0.10, 0.86, 0.85, 0.98, "NDC");
    pt_pT->SetFillColor(kWhite);
    pt_pT->SetBorderSize(0);
    pt_pT->SetTextAlign(12);
    pt_pT->AddText(Form("Transverse Momentum (p_T) - Event %d with R= %.1f and ptmin= %.1f GeV Jets", eventToShow, R, ptmin));
    pt_pT->AddText(Form("Underflow (#eta < -5): %.2f GeV  |  Overflow (#eta > 5): %.2f GeV", pt_under, pt_over));
    pt_pT->Draw();

    for (unsigned i = 0; i < jets.size(); i++) {
        TEllipse *circle = new TEllipse(jets[i].rap(), jets[i].phi_std(), R, R);
        circle->SetFillStyle(0); 
        circle->SetLineColor(kRed);
        circle->SetLineWidth(3);
        circle->Draw("same");
    }

    c1->SaveAs(Form("Event_%d_Jets.png", eventToShow));
    c1->SaveAs(Form("Event_%d_Jets.pdf", eventToShow));
    
    cout << "Plots saved as Event_" << eventToShow << "_Jets.png and .pdf" << endl;
    
    f->Close();
    return 0;
}