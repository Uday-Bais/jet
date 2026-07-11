#include <iostream>
#include "Pythia8/Pythia.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;

int main() {
    int nevents = 1000000;
    Pythia8::Pythia pythia;
    
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000");
    pythia.readString("HardQCD:all = on");
    
    // added this later as no pi zero was found

    //pythia.readString("111:mayDecay = off");
    
    pythia.init();


    TFile *outputFile = new TFile("pythia_events_1m_pt.root", "RECREATE");
    TTree *tree = new TTree("ParticleTree", "Kinematics of the produced particles");
    
    int EventID, ParticleID;
    double pT, eta, phi, mass, Px, Py, Pz, Energy;
    
    tree->Branch("EventID",    &EventID,    "EventID/I");
    tree->Branch("ParticleID", &ParticleID, "ParticleID/I");
    tree->Branch("pT",         &pT,         "pT/D");
    tree->Branch("eta",        &eta,        "eta/D");
    tree->Branch("phi",        &phi,        "phi/D");
    tree->Branch("mass",       &mass,       "mass/D");
    tree->Branch("Px",         &Px,         "Px/D");
    tree->Branch("Py",         &Py,         "Py/D");
    tree->Branch("Pz",         &Pz,         "Pz/D");
    tree->Branch("Energy",     &Energy,     "Energy/D");
    
    for(int i = 0; i < nevents; i++) {
        if(!pythia.next()) continue;
        int entries = pythia.event.size();
        
        for(int j = 0; j < entries; j++) {

            if (!pythia.event[j].isFinal()) continue;
            if (!pythia.event[j].isVisible()) continue;
            if (pythia.event[j].pT() <= 0.5) continue;
            EventID    = i;
            ParticleID = pythia.event[j].id();
            pT         = pythia.event[j].pT();
            eta        = pythia.event[j].eta();
            phi        = pythia.event[j].phi();
            mass       = pythia.event[j].m();
            Px         = pythia.event[j].px();
            Py         = pythia.event[j].py();
            Pz         = pythia.event[j].pz();
            Energy     = pythia.event[j].e();
            

            tree->Fill();
        }
    }
    

    tree->Write();
    outputFile->Close();
    
    cout << "successfully saved" << endl;
    return 0;
}
