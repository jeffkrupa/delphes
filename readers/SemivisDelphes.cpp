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

static int NMAX = 20;

//---------------------------------------------------------------------------

struct FJ
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float m = 0;
};

//---------------------------------------------------------------------------

template <typename T>
void
fill(vector<float> &vattr, vector<PFCand> &particles, T fn_attr)
{
  vattr.clear();
  for (auto& p : particles)
    vattr.push_back(fn_attr(p));
}


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
  TBranch* fjbranch = (TBranch*)itree->GetBranch("FatJet");
  std::cout << "NEVT: " << nevt << std::endl;

  vector<float> fjpt, fjeta, fjphi, fje;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX);

  float pfmet=-99., pfmetphi=-99., puppimet=-99., puppimetphi=-99.;
  float ht=-99., float fjDEta=-99.;

  // branches
  TBranch* b_pfmet = tout->Branch("pfmet",&pfmet, "pfmet/F");
  TBranch* b_pfmetphi = tout->Branch("pfmetphi",&pfmetphi, "pfmetphi/F");
  TBranch* b_puppimet = tout->Branch("puppimet",&puppimet, "puppimet/F");
  TBranch* b_puppimetphi = tout->Branch("puppimetphi",&puppimetphi, "puppimetphi/F");
  TBranch* b_ht = tout->Branch("ht",&ht, "ht/F");
  TBranch* b_fjDEta = tout->Branch("fjDEta",&fjDEta, "fjDEta/F");
  
  // fat jet branches
  tout->Branch("fjPt", &fjpt);
  tout->Branch("fjEta", &fjeta);
  tout->Branch("fjPhi", &fjphi);
  tout->Branch("fjM", &fje);

  ExRootProgressBar progressBar(nevt);
  
  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);

    if (k%100==0)
      std::cout << k << " / " << nevt << std::endl;

    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);

    for (unsigned int j=0; j<npfs; j++){
      FJ tmpfj;
      tmpfj.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      tmpfj.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmpfj.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmpfj.m = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.puppi = itree->GetLeaf("ParticleFlowCandidate.PuppiW")->GetValue(j);
      tmppf.hardfrac = itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j);
      if (itree->GetLeaf("ParticleFlowCandidate.Charge")->GetValue(j)!=0){
        if (itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j)==1)
          tmppf.vtxid = 0;
        else
          tmppf.vtxid = 1;
      }
      else
        tmppf.vtxid = -1;
      input_particles.push_back(tmppf);
    }

    tout->Fill();
    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
