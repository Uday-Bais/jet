#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> //for std::sort

// ROOT headers
#include "TFile.h"
#include "TTree.h"

using namespace std;

// Struct to represent a particle/pseudojet
struct Particle {
    int id;
    double pt;
    double eta;
    double phi;
};

// Calculate delta R squared
double deltaR2(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    
    while (dphi > M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    
    return deta * deta + dphi * dphi;
}

int main() {
    TFile *inputFile = new TFile("pythia_events_1m_pt.root", "READ");
    if (!inputFile->IsOpen()) {
        cerr << "Error: Could not open ROOT file!" << endl;
        return 1;
    }
    
    TTree *tree = (TTree*)inputFile->Get("ParticleTree");
    
    int EventID;
    double pT, eta, phi;
    
    tree->SetBranchAddress("EventID", &EventID);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("phi", &phi);
    
    int target_event = 10;
    vector<Particle> particles;
    
    // Vector to store the final jets
    vector<Particle> final_jets; 
    
    int original_id_counter = 0;
    
    // pT cut 
    double pt_cut = 0.5; 
    
    int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (EventID > target_event) break; 
        
        if (EventID == target_event) {
            if (pT <= pt_cut) continue; 
            
            Particle p;
            p.id = original_id_counter++;
            p.pt = pT;
            p.eta = eta;
            p.phi = phi;
            particles.push_back(p);
        }
    }
    inputFile->Close();
    
    if (particles.empty()) {
        cout << "No particles left " << target_event << " with pT > " << pt_cut << " GeV." << endl;
        return 1;
    }
    
    cout << "Loaded " << particles.size() << " particles for Event " << target_event << " (pT > " << pt_cut << " GeV)\n";
    cout << "===========================================================\n";

    double R = 1;
    double R2 = R * R;
    int step = 1;
    int next_id = 10000; 
    
    cout << fixed << setprecision(6);
    
    while (!particles.empty()) {
        cout << "\n====================== STEP " << step << " ======================\n";
        cout << "Active Particles: " << particles.size() << "\n\n";
        
        double min_dist = 1e15; 
        int min_i = -1;
        int min_j = -1;
        bool is_beam = false;
        
        // 1. Calculate and print all d_iB
        cout << "--- Particle-to-Beam Distances (d_iB) ---\n";
        for (size_t i = 0; i < particles.size(); i++) {
            double inv_pt2 = 1.0 / (particles[i].pt * particles[i].pt);
            double d_iB = inv_pt2;
            
            cout << "d_{" << particles[i].id << "B} = " << d_iB << "\n";
            
            if (d_iB < min_dist) {
                min_dist = d_iB;
                min_i = i;
                is_beam = true;
            }
        }
        
        cout << "\n--- Particle-to-Particle Distances (d_ij) ---\n";
        // 2. Calculate and print all d_ij
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i + 1; j < particles.size(); j++) {
                double inv_pti2 = 1.0 / (particles[i].pt * particles[i].pt);
                double inv_ptj2 = 1.0 / (particles[j].pt * particles[j].pt);
                
                double pt_min2 = min(inv_pti2, inv_ptj2);
                double dR2 = deltaR2(particles[i].eta, particles[i].phi, particles[j].eta, particles[j].phi);
                
                double d_ij = pt_min2 * (dR2 / R2);
                
                cout << "d_{" << particles[i].id << "," << particles[j].id << "} = " << d_ij 
                     << "  [min(1/pT^2) = " << pt_min2 << ", dR^2/R^2 = " << (dR2/R2) << "]\n";
                
                if (d_ij < min_dist) {
                    min_dist = d_ij;
                    min_i = i;
                    min_j = j;
                    is_beam = false;
                }
            }
        }
        
        // 3. Print the decision 
        cout << "\n--- STEP " << step << " ACTION DECISION ---\n";
        
        if (is_beam) {
            cout << "ABSOLUTE MINIMUM: d_{" << particles[min_i].id << "B} = " << min_dist << "\n";
            cout << "ACTION: Particle " << particles[min_i].id << " is isolated. Declaring FINAL JET.\n";
            cout << "    Jet Kinematics: pT = " << particles[min_i].pt 
                 << ", eta = " << particles[min_i].eta 
                 << ", phi = " << particles[min_i].phi << "\n";
            
            
            final_jets.push_back(particles[min_i]);
            
            particles.erase(particles.begin() + min_i);
        } else {
            cout << "ABSOLUTE MINIMUM: d_{" << particles[min_i].id << "," << particles[min_j].id << "} = " << min_dist << "\n";
            cout << "ACTION: Merging Particle " << particles[min_i].id << " and Particle " << particles[min_j].id << " into New Pseudojet " << next_id << ".\n";
            
            Particle merged;
            merged.id = next_id++;
            merged.pt = particles[min_i].pt + particles[min_j].pt;
            merged.eta = (particles[min_i].pt * particles[min_i].eta + particles[min_j].pt * particles[min_j].eta) / merged.pt;
            merged.phi = (particles[min_i].pt * particles[min_i].phi + particles[min_j].pt * particles[min_j].phi) / merged.pt;
            
            cout << "    final properties: pT = " << merged.pt 
                 << ", eta = " << merged.eta 
                 << ", phi = " << merged.phi << "\n";
            
            particles.erase(particles.begin() + max(min_i, min_j));
            particles.erase(particles.begin() + min(min_i, min_j));
            
            particles.push_back(merged);
        }
        
        step++;
    }
    
    // Final Summary
    cout << "\n===========================================================\n";
    cout << "COMPLETE. No active particles remaining.\n";
    cout << "===========================================================\n\n";
    
    // Sort final jets by pT descending
    sort(final_jets.begin(), final_jets.end(), [](const Particle& a, const Particle& b) {
        return a.pt > b.pt;
    });

    cout << "FINAL JET SUMMARY\n";
    cout << "Total Final Jets Extracted: " << final_jets.size() << "\n\n";
    cout << left << setw(10) << "Jet ID" 
         << setw(15) << "pT [GeV]" 
         << setw(15) << "eta" 
         << setw(15) << "phi" << "\n";
    cout << "------------------------------------------------------\n";
    
    for (const auto& jet : final_jets) {
        cout << left << setw(10) << jet.id 
             << setw(15) << jet.pt 
             << setw(15) << jet.eta 
             << setw(15) << jet.phi << "\n";
    }
    
    return 0;
}
