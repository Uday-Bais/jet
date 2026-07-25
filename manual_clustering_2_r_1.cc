#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> 
#include <fstream>

#include "TFile.h"
#include "TTree.h"

using namespace std;

struct Particle {
    int id;
    double pt;
    double y;      
    double phi;
    double px;
    double py;
    double pz;
    double e;
};

double deltaR2(double y1, double phi1, double y2, double phi2) {
    double dy = y1 - y2;
    double dphi = phi1 - phi2;
    
    while (dphi > M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    
    return dy * dy + dphi * dphi;
}

double safe_rapidity(double e, double pz) {
    double e_plus_pz = e + pz;
    double e_minus_pz = e - pz;
    
    if (e_minus_pz <= 1e-9) {
        return 10.0; 
    } else if (e_plus_pz <= 1e-9) {
        return -10.0; 
    } else {
        return 0.5 * log(e_plus_pz / e_minus_pz);
    }
}

int main() {
    TFile *inputFile = new TFile("full_history_pythia_20gev_ptHat.root", "READ");
    if (!inputFile->IsOpen()) {
        cerr << "Error: Could not open ROOT file!" << endl;
        return 1;
    }
    
    ofstream outFile("jet_output.txt");
    if (!outFile.is_open()) {
        cerr << "Error: Could not create output text file!" << endl;
        return 1;
    }
    
    TTree *tree = (TTree*)inputFile->Get("ParticleTree");
    
    int EventID = 0;
    double pT = 0.0, phi = 0.0;
    double Energy = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0;
    bool IsFinal = false;
    int pdgID = 0;
    
    tree->SetBranchAddress("EventID", &EventID);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("Energy", &Energy);
    tree->SetBranchAddress("Px", &Px);
    tree->SetBranchAddress("Py", &Py);
    tree->SetBranchAddress("Pz", &Pz);
    tree->SetBranchAddress("IsFinal", &IsFinal);
    tree->SetBranchAddress("pdgID", &pdgID);
    
    int target_event = 0;
    vector<Particle> particles;
    vector<Particle> final_jets; 
    
    int original_id_counter = 0;
    
    int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (EventID > target_event) break; 
        
        if (EventID == target_event) {
            if (!IsFinal) continue;
            
            int abs_pdg = std::abs(pdgID);
            if (abs_pdg == 12 || abs_pdg == 14 || abs_pdg == 16) continue;
            
            Particle p;
            p.id = original_id_counter++;
            p.pt = pT;
            p.phi = phi;
            p.px = Px;
            p.py = Py;
            p.pz = Pz;
            p.e  = Energy; 
            
            p.y = safe_rapidity(p.e, p.pz);
            
            particles.push_back(p);
        }
    }
    inputFile->Close();
    
    if (particles.empty()) {
        outFile << "No visible final-state particles loaded for Event " << target_event << "\n";
        outFile.close();
        return 1;
    }
    
    outFile << "Loaded " << particles.size() << " particles for Event " << target_event << "\n";
    outFile << "===========================================================\n";

    double R = 1;
    double R2 = R * R;
    int step = 1;
    int next_id = 10000; 
    
    outFile << fixed << setprecision(6);
    
    while (!particles.empty()) {
        double min_dist = -1.0; 
        int min_i = -1;
        int min_j = -1;
        bool is_beam = false;
        bool first_calc = true; 
        
        for (size_t i = 0; i < particles.size(); i++) {
            double pt = particles[i].pt;
            double inv_pt2 = (pt > 1e-10) ? 1.0 / (pt * pt) : 1e20; 
            
            if (first_calc || inv_pt2 < min_dist) {
                min_dist = inv_pt2;
                min_i = i;
                is_beam = true;
                first_calc = false;
            }
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i + 1; j < particles.size(); j++) {
                double pti = particles[i].pt;
                double ptj = particles[j].pt;
                
                double inv_pti2 = (pti > 1e-10) ? 1.0 / (pti * pti) : 1e20;
                double inv_ptj2 = (ptj > 1e-10) ? 1.0 / (ptj * ptj) : 1e20;
                
                double pt_min2 = min(inv_pti2, inv_ptj2);
                
                double dR2 = deltaR2(particles[i].y, particles[i].phi, particles[j].y, particles[j].phi);
                
                double d_ij = pt_min2 * (dR2 / R2);
                
                if (first_calc || d_ij < min_dist) {
                    min_dist = d_ij;
                    min_i = i;
                    min_j = j;
                    is_beam = false;
                    first_calc = false;
                }
            }
        }
        
        if (is_beam) {
            final_jets.push_back(particles[min_i]);
            particles.erase(particles.begin() + min_i);
        } else {
            Particle merged;
            merged.id = next_id++;
            
            merged.px = particles[min_i].px + particles[min_j].px;
            merged.py = particles[min_i].py + particles[min_j].py;
            merged.pz = particles[min_i].pz + particles[min_j].pz;
            merged.e  = particles[min_i].e + particles[min_j].e;
            
            merged.pt = sqrt(merged.px * merged.px + merged.py * merged.py);
            merged.phi = atan2(merged.py, merged.px);
            
            merged.y = safe_rapidity(merged.e, merged.pz);
            
            particles.erase(particles.begin() + max(min_i, min_j));
            particles.erase(particles.begin() + min(min_i, min_j));
            
            particles.push_back(merged);
        }
        step++;
    }
    
    sort(final_jets.begin(), final_jets.end(), [](const Particle& a, const Particle& b) {
        return a.pt > b.pt;
    });

    outFile << "\nFINAL JET SUMMARY (Filtered for pT > 10 GeV)\n";
    outFile << left << setw(10) << "Jet ID" 
            << setw(12) << "pT" 
            << setw(12) << "y" 
            << setw(12) << "phi" 
            << setw(12) << "pX"
            << setw(12) << "pY"
            << setw(12) << "pZ"
            << setw(12) << "E" << "\n";
    outFile << "--------------------------------------------------------------------------------------------\n";
    
    double jet_pt_cut = 10.0;
    int valid_jets = 0;
    
    for (const auto& jet : final_jets) {
        if (jet.pt > jet_pt_cut) {
            outFile << left << setw(10) << jet.id 
                    << setw(12) << jet.pt 
                    << setw(12) << jet.y 
                    << setw(12) << jet.phi 
                    << setw(12) << jet.px
                    << setw(12) << jet.py
                    << setw(12) << jet.pz
                    << setw(12) << jet.e << "\n";
            valid_jets++;
        }
    }
    
    outFile << "\nTotal Jets above threshold: " << valid_jets << "\n";
    
    outFile.close();
    
    return 0;
}
