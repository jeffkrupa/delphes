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


//---------------------------------------------------------------------------


struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float puppi = 1;
  float hardfrac = 1;  
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

  float genjetpt=-99., genjeteta=-99., genjetphi=-99., genjetm=-99.;
  float puppijetpt=-99., puppijeteta=-99., puppijetphi=-99., puppijetm=-99.;
  float truthjetpt=-99., truthjeteta=-99., truthjetphi=-99., truthjetm=-99.;
  
  TBranch* b_genjetpt = tout->Branch("genjetpt",&genjetpt, "genjetpt/F");
  TBranch* b_genjeteta = tout->Branch("genjeteta",&genjeteta, "genjeteta/F");
  TBranch* b_genjetphi = tout->Branch("genjetphi",&genjetphi, "genjetphi/F");
  TBranch* b_genjetm = tout->Branch("genjetm",&genjetm, "genjetm/F");
  TBranch* b_puppijetpt = tout->Branch("puppijetpt",&puppijetpt, "puppijetpt/F");
  TBranch* b_puppijeteta = tout->Branch("puppijeteta",&puppijeteta, "puppijeteta/F");
  TBranch* b_puppijetphi = tout->Branch("puppijetphi",&puppijetphi, "puppijetphi/F");
  TBranch* b_puppijetm = tout->Branch("puppijetm",&puppijetm, "puppijetm/F");
  TBranch* b_truthjetpt = tout->Branch("truthjetpt",&truthjetpt, "truthjetpt/F");
  TBranch* b_truthjeteta = tout->Branch("truthjeteta",&truthjeteta, "truthjeteta/F");
  TBranch* b_truthjetphi = tout->Branch("truthjetphi",&truthjetphi, "truthjetphi/F");
  TBranch* b_truthjetm = tout->Branch("truthjetm",&truthjetm, "truthjetm/F");

  ExRootProgressBar progressBar(nevt);
  
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  //int activeAreaRepeats = 1;
  //double ghostArea = 0.01;
  //double ghostEtaMax = 7.0;

  //fastjet::GhostedAreaSpec *activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  //fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);
  fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.8);

  double radius = 0.8;
  double sdZcut = 0.1;
  double sdBeta = 0.;
  
  fastjet::contrib::SoftDrop softDrop = fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);

  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);
    //if (k>100)
    //break;

    //if (k%100==0)
    //std::cout << k << " / " << nevt << std::endl;

    TLorentzVector higgs;
    higgs.SetPtEtaPhiE(0.,0.,0.,0.);
    unsigned int nparts = partbranch->GetEntries();
    nparts = itree->GetLeaf("Particle_size")->GetValue(0);
    for (unsigned int j=0; j<nparts; j++){
      if (itree->GetLeaf("Particle.PID")->GetValue(j)==25){
	higgs.SetPtEtaPhiE(itree->GetLeaf("Particle.PT")->GetValue(j),itree->GetLeaf("Particle.Eta")->GetValue(j),itree->GetLeaf("Particle.Phi")->GetValue(j),itree->GetLeaf("Particle.E")->GetValue(j));
	break;
      }
    }
    
    unsigned int ngenjets = genjetbranch->GetEntries();
    ngenjets = itree->GetLeaf("GenJet_size")->GetValue(0);
    for (unsigned int j=0; j<ngenjets; j++){
      TLorentzVector tmpjet;
      tmpjet.SetPtEtaPhiM(itree->GetLeaf("GenJet.PT")->GetValue(j),itree->GetLeaf("GenJet.Eta")->GetValue(j),itree->GetLeaf("GenJet.Phi")->GetValue(j),itree->GetLeaf("GenJet.Mass")->GetValue(j));
      if (tmpjet.DeltaR(higgs)<0.8){
	genjetpt = tmpjet.Pt();
	genjeteta = tmpjet.Eta();
	genjetphi = tmpjet.Phi();
	genjetm = tmpjet.M();
	break;
      }
    }

    input_particles.clear();
    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);

    for (unsigned int j=0; j<npfs; j++){
      PFCand tmppf;
      tmppf.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      tmppf.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmppf.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmppf.e = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.puppi = itree->GetLeaf("ParticleFlowCandidate.PuppiW")->GetValue(j);
      tmppf.hardfrac = itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j);
      input_particles.push_back(tmppf);
    }

    // sorting input particles by pT
    sort(input_particles.begin(), input_particles.end(), comp_p4);

    vector<fastjet::PseudoJet> finalStates_puppi;
    vector<fastjet::PseudoJet> finalStates_truth;
    vector<fastjet::PseudoJet> finalStates_pf;
    for(auto &p : input_particles){
      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(p.pt,p.eta,p.phi,p.e);
      finalStates_pf.emplace_back(tmp.Px(), tmp.Py(), tmp.Pz(), tmp.E());
      if (p.puppi>0.)
	finalStates_puppi.emplace_back(p.puppi*tmp.Px(), p.puppi*tmp.Py(), p.puppi*tmp.Pz(), p.puppi*tmp.E());
      if (p.hardfrac>0.)
	finalStates_truth.emplace_back(p.hardfrac*tmp.Px(), p.hardfrac*tmp.Py(), p.hardfrac*tmp.Pz(), p.hardfrac*tmp.E());
    }

    fastjet::ClusterSequence seq_puppi(sorted_by_pt(finalStates_puppi), *jetDef);
    fastjet::ClusterSequence seq_truth(sorted_by_pt(finalStates_truth), *jetDef);
    //fastjet::ClusterSequence seq_pf(sorted_by_pt(finalStates_pf), *jetDef);

    vector<fastjet::PseudoJet> allJets_puppi(sorted_by_pt(seq_puppi.inclusive_jets(400.)));
    vector<fastjet::PseudoJet> allJets_truth(sorted_by_pt(seq_truth.inclusive_jets(400.)));
    //vector<fastjet::PseudoJet> allJets_pf(sorted_by_pt(seq_pf.inclusive_jets(400.)));

    //for (auto& pfJet : allJets_pf) {
    //if (pfJet.m()<0 && pfJet.perp() > 400.)
    //	std::cout << "WTF: " << pfJet.m() << std::endl;
    //}

    for (auto& puppiJet : allJets_puppi) {
      if (puppiJet.perp() < 400.)
	break;

      fastjet::PseudoJet sdJet = (softDrop)(puppiJet);

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(puppiJet.perp(),puppiJet.eta(),puppiJet.phi(),sdJet.m());
      if (tmp.DeltaR(higgs)<0.8){
	puppijetpt = tmp.Pt();
	puppijeteta = tmp.Eta();
	puppijetphi = tmp.Phi();
	puppijetm = tmp.M();
	break;
      }
    }

    for (auto& truthJet : allJets_truth) {
      if (truthJet.perp() < 400.)
	break;

      fastjet::PseudoJet sdJet = (softDrop)(truthJet);

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(truthJet.perp(),truthJet.eta(),truthJet.phi(),sdJet.m());
      if (tmp.DeltaR(higgs)<0.8){
	truthjetpt = tmp.Pt();
	truthjeteta = tmp.Eta();
	truthjetphi = tmp.Phi();
	truthjetm = tmp.M();
	break;
      }
    }

    tout->Fill();

    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
