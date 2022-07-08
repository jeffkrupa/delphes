#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <array>
#include <stdlib.h>
#include <functional>
#include <time.h>
#include <math.h>
#include <numeric>

#include "TROOT.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

static int NMAX = 1000;

//---------------------------------------------------------------------------


struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float pdgid = 0;
  float d0 = 0;
  float dz = 0;
  float charge = 0;
};



//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  srand(time(NULL));

  if(argc < 3) {
    cout << " Usage: " << "PapuDelphes" << " input_file"
         << " output_file" << endl;
    cout << " input_file - input file in ROOT format," << endl;
    cout << " output_file - output file in ROOT format" << endl;
    return 1;
  }

  // figure out how to read the file here 
  //

  TFile* ifile = TFile::Open(argv[1], "READ");
  TTree* itree = (TTree*)ifile->Get("Delphes;1");

  auto* fout = TFile::Open(argv[2], "RECREATE");
  auto* tout = new TTree("events", "events");

  unsigned int nevt = itree->GetEntries();
  TBranch* partbranch = (TBranch*)itree->GetBranch("Particle");
  TBranch* pfbranch = (TBranch*)itree->GetBranch("ParticleFlowCandidate");
  TBranch* genjetbranch = (TBranch*)itree->GetBranch("GenJet");
  std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;

  vector<PFCand> output_particles;
  output_particles.reserve(NMAX);

  vector<float> vpt, veta, vphi, ve, vpuppi, vpdgid, vhardfrac, vcharge, vd0, vdz;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX); vpdgid.reserve(NMAX); 
  vcharge.reserve(NMAX);
  vd0.reserve(NMAX);
  vdz.reserve(NMAX);

  float genjetpt=-99., genjeteta=-99., genjetphi=-99., genjetm=-99.;
  float puppijetpt=-99., puppijeteta=-99., puppijetphi=-99., puppijetm=-99.;
  float truthjetpt=-99., truthjeteta=-99., truthjetphi=-99., truthjetm=-99.;
  float pfjetpt=-99., pfjeteta=-99., pfjetphi=-99., pfjetm=-99.;

  // jet branches
  // PF cand branches
  tout->Branch("pt", &vpt);
  tout->Branch("eta", &veta);
  tout->Branch("phi", &vphi);
  tout->Branch("e", &ve);
  tout->Branch("charge", &vcharge);
  tout->Branch("d0", &vd0);
  tout->Branch("dz", &vdz);
  tout->Branch("pdgid", &vpdgid);

  ExRootProgressBar progressBar(nevt);
  
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.8);

  // Soft drop -- only needed for large radius jets
  double radius = 0.8;
  double sdZcut = 0.1;
  double sdBeta = 0.;
  
  fastjet::contrib::SoftDrop softDrop = fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);

  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);

    if (k%100==0)
      std::cout << k << " / " << nevt << std::endl;

    input_particles.clear();
    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);

    for (unsigned int j=0; j<npfs; j++){
      PFCand tmppf;
      tmppf.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      tmppf.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmppf.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmppf.e = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.pdgid = itree->GetLeaf("ParticleFlowCandidate.PID")->GetValue(j);
      tmppf.charge = itree->GetLeaf("ParticleFlowCandidate.Charge")->GetValue(j);
      tmppf.d0 = 0.;
      if (abs(itree->GetLeaf("ParticleFlowCandidate.D0")->GetValue(j) > 1e-5)){
	tmppf.d0 = itree->GetLeaf("ParticleFlowCandidate.D0")->GetValue(j);
      }	
      else{
	tmppf.d0 = 0;
      }
      tmppf.dz = itree->GetLeaf("ParticleFlowCandidate.DZ")->GetValue(j);

      input_particles.push_back(tmppf);
    }

    // sorting input particles by pT
    sort(input_particles.begin(), input_particles.end(), comp_p4);

    vector<fastjet::PseudoJet> finalStates;
    int pfid = 0;
    for(auto &p : input_particles){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(p.pt,p.eta,p.phi,p.e);

      fastjet::PseudoJet curjet(tmp.Px(), tmp.Py(), tmp.Pz(), tmp.E());
      curjet.set_user_index(pfid);    
      finalStates.emplace_back(curjet);
      
      pfid += 1;
    }

    fastjet::ClusterSequence seq(sorted_by_pt(finalStates), *jetDef);

    float minpt = 400;
    vector<fastjet::PseudoJet> allJets(sorted_by_pt(seq.inclusive_jets(minpt)));

    bool has_match = false;

    /*

    for (auto& jet : allJets) {
      if (jet.perp() < minpt)
	break;

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(jet.perp(),jet.eta(),jet.phi(),jet.m());
      if ((tmp.DeltaR(Higgs)<0.8) && (tmp.DeltaR(b1)<0.8) && (tmp.DeltaR(b2)<0.8)){
	has_match = true;

      }
    }

    */

    tout->Fill();

    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
