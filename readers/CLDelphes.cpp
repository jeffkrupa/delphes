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
#include <fstream>
#include "TROOT.h"
#include <cassert>
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
//#include "fastjet/EnergyCorrelator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/contribs/EnergyCorrelator/EnergyCorrelator.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

static int NMAX = 200;

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


template <typename T>
void
fill(vector<float> &vattr, vector<PFCand> &particles, T fn_attr)
{
  vattr.clear();
  for (auto& p : particles)
    vattr.push_back(fn_attr(p));
}

//---------------------------------------------------------------------------
//
Double_t deltaPhi(float phi1, float phi2){
    Double_t deltaPhi = TMath::Abs(phi1-phi2);
    if(deltaPhi > TMath::Pi())
        deltaPhi = TMath::TwoPi() - deltaPhi;
    return deltaPhi;
}

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
  //std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;
  vector<PFCand> output_particles;
  output_particles.reserve(NMAX);

  vector<float> vpt, veta, vphi, ve, vpuppi, vpdgid, vhardfrac, vcharge, vd0, vdz;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX); vpdgid.reserve(NMAX); 
  vcharge.reserve(NMAX);
  vd0.reserve(NMAX);
  vdz.reserve(NMAX);

  float jettype=-1.; // 0: gluon, 1: light quark, 2: charm, 3: b, 4: Higgs, 5: Z, 6: top
  float m_parton_pt = 0.;
  float m_parton_eta = 0.;
  float m_parton_phi = 0.;
  float m_parton_e = 0.;
  float dau1_parton_pt = 0.;
  float dau1_parton_eta = 0.;
  float dau1_parton_phi = 0.;
  float dau1_parton_e = 0.;
  float dau2_parton_pt = 0.;
  float dau2_parton_eta = 0.;
  float dau2_parton_phi = 0.;
  float dau2_parton_e = 0.;
  float jet_pt = 0.;
  float jet_eta = 0.;
  float jet_phi = 0.;
  float jet_e = 0.;
  float jet_msd = 0.;
  float jet_n2 = -99.;

  TBranch* b_m_parton_pt = tout->Branch("m_parton_pt",&m_parton_pt, "m_parton_pt/F");
  TBranch* b_m_parton_eta = tout->Branch("m_parton_eta",&m_parton_eta, "m_parton_eta/F");
  TBranch* b_m_parton_phi = tout->Branch("m_parton_phi",&m_parton_phi, "m_parton_phi/F");
  TBranch* b_m_parton_e = tout->Branch("m_parton_e",&m_parton_e, "m_parton_e/F");

  TBranch* b_jettype = tout->Branch("jettype",&jettype, "jettype/F");
  TBranch* b_jet_pt = tout->Branch("jet_pt",&jet_pt, "jet_pt/F");
  TBranch* b_jet_eta = tout->Branch("jet_eta",&jet_eta, "jet_eta/F");
  TBranch* b_jet_phi = tout->Branch("jet_phi",&jet_phi, "jet_phi/F");
  TBranch* b_jet_e = tout->Branch("jet_e",&jet_e, "jet_e/F");
  TBranch* b_jet_msd = tout->Branch("jet_msd",&jet_msd, "jet_msd/F");
  TBranch* b_jet_n2 = tout->Branch("jet_n2",&jet_n2, "jet_n2/F");

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
  EnergyCorrelatorNseries N2fj(2,1.0,EnergyCorrelator::pt_R);

  //std::ofstream outfile;
  //outfile.open("../../jettype.txt", std::ios_base::trunc);

  for (unsigned int k=0; k<100; k++){
    //std::cout <<"NEW EVENT ========" << std::endl;
    itree->GetEntry(k);
    //if( k>10) break;
    if (k%1000==0)
      //std::cout << k << " / " << nevt << std::endl;

    input_particles.clear();
    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);
    std::cout << npfs << std::endl;
    for (unsigned int j=0; j<npfs; j++){
      PFCand tmppf;
      tmppf.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      if (tmppf.pt < 1.0) { continue; }
      tmppf.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmppf.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmppf.e = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.pdgid = itree->GetLeaf("ParticleFlowCandidate.PID")->GetValue(j);
      tmppf.charge = itree->GetLeaf("ParticleFlowCandidate.Charge")->GetValue(j);
      tmppf.d0 = 0.;
      if (abs(itree->GetLeaf("ParticleFlowCandidate.D0")->GetValue(j) > 1e-6)){
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

    float minpt = 500;
    float maxeta = 2.5;
    vector<fastjet::PseudoJet> allJets(sorted_by_pt(seq.inclusive_jets(minpt)));
    std::cout << "allJets.size()" << allJets.size() << std::endl;
    bool has_match = false;
    int jetidx = 0;
    for (auto& jet : allJets) {
      
      //std::cout << "JET " << jetidx << " PT:"<< jet.perp() << " ETA:" << jet.eta() << std::endl;
      ++jetidx;
      if (jet.perp() < minpt || abs(jet.eta()) > maxeta )
	break;

      output_particles.clear();


      jettype = -1.;
      m_parton_pt = 0.;
      m_parton_eta = 0.;
      m_parton_phi = 0.;
      m_parton_e = 0.;
      dau1_parton_pt = 0.;
      dau1_parton_eta = 0.;
      dau1_parton_phi = 0.;
      dau1_parton_e = 0.;
      dau2_parton_pt = 0.;
      dau2_parton_eta = 0.;
      dau2_parton_phi = 0.;
      dau2_parton_e = 0.;

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(jet.perp(),jet.eta(),jet.phi(),jet.m());

      unsigned int nparts = partbranch->GetEntries();
      nparts = itree->GetLeaf("Particle_size")->GetValue(0);

      bool has_zprime = false;
      TLorentzVector zprime(0.,0.,0.,0);
      //std::cout <<"here0"<<std::endl;
      for (unsigned int w=0; w<nparts; w++){
	//std::cout << "mass:" << itree->GetLeaf("Particle.Mass")->GetValue(w) << std::endl;
	if (itree->GetLeaf("Particle.PID")->GetValue(w) == 55){
	  has_zprime = true;
	  zprime.SetPtEtaPhiM(itree->GetLeaf("Particle.PT")->GetValue(w),itree->GetLeaf("Particle.Eta")->GetValue(w),itree->GetLeaf("Particle.Phi")->GetValue(w),itree->GetLeaf("Particle.Mass")->GetValue(w));
	  m_parton_pt = itree->GetLeaf("Particle.PT")->GetValue(w);
	  m_parton_eta = itree->GetLeaf("Particle.Eta")->GetValue(w);
	  m_parton_phi = itree->GetLeaf("Particle.Phi")->GetValue(w);
	  m_parton_e = itree->GetLeaf("Particle.E")->GetValue(w);
	  break;
	}
      }
      /*
      std::cout << "start_particles" << std::endl;
      for (unsigned int w=0; w<nparts; w++){
          std::cout << "particle" << w << " pdgId=" << itree->GetLeaf("Particle.PID")->GetValue(w) << " M1="<<itree->GetLeaf("Particle.M1")->GetValue(w) << " M2=" << itree->GetLeaf("Particle.M2")->GetValue(w) << " D1=" << itree->GetLeaf("Particle.D1")->GetValue(w) << " D2=" << itree->GetLeaf("Particle.D2")->GetValue(w) << std::endl;
      }
      std::cout << "end_particles" << std::endl;
      std::cout << has_zprime << std::endl;
      //std::cout <<"here1"<<std::endl;
      */
      if (has_zprime){	
	TLorentzVector q1(0.,0.,0.,0);
	TLorentzVector q2(0.,0.,0.,0);

	int qsfound = 0;
        int q1idx = nparts-1;
        int q2idx = nparts-2;
        q1.SetPtEtaPhiE(itree->GetLeaf("Particle.PT")->GetValue(q1idx),itree->GetLeaf("Particle.Eta")->GetValue(q1idx),itree->GetLeaf("Particle.Phi")->GetValue(q1idx),itree->GetLeaf("Particle.E")->GetValue(q1idx));
        q2.SetPtEtaPhiE(itree->GetLeaf("Particle.PT")->GetValue(q2idx),itree->GetLeaf("Particle.Eta")->GetValue(q2idx),itree->GetLeaf("Particle.Phi")->GetValue(q2idx),itree->GetLeaf("Particle.E")->GetValue(q2idx));
        assert(itree->GetLeaf("Particle.M1")->GetValue(q1idx) == itree->GetLeaf("Particle.M1")->GetValue(q2idx) );
        //std::cout << "\tq1 : pt=" << q1.Pt()  << " eta=" << q1.Eta() << " phi=" << q1.Phi() << " pgdId=" << itree->GetLeaf("Particle.PID")->GetValue(q1idx) << std::endl;
        //std::cout << "\tq2 : pt=" << q2.Pt()  << " eta=" << q2.Eta() << " phi=" << q2.Phi() << " pgdId=" << itree->GetLeaf("Particle.PID")->GetValue(q2idx) << std::endl;
        //std::cout << tmp.DeltaR(zprime) << "/" << tmp.DeltaR(q1) << "/" << tmp.DeltaR(q2) << std::endl;
	if ((tmp.DeltaR(zprime)<0.5) && (tmp.DeltaR(q1)<0.8) && (tmp.DeltaR(q2)<0.8)) {
	  jettype = 1.;
	}
      }
      else {
	  jettype = 0.; 
      } 
      
      //std::cout <<"here2"<<std::endl;
      if (jettype>-1.){
	//std::cout << jettype << std::endl;
        //outfile << jettype << "\n";

	fastjet::PseudoJet sdJet = (softDrop)(jet);

	jet_pt = tmp.Pt();
	jet_eta = tmp.Eta();
	jet_phi = tmp.Phi();
	jet_e = tmp.E();
	jet_msd = sdJet.m();

	//cout << N2fj(jet) << endl;

	jet_n2 = N2fj(jet);
	//cout << "Before/after softdrop: " << jet.constituents().size() << " / " << sdJet.constituents().size() << endl;

	// fill constituents
        //std::cout <<jet.constituents().size() << std::endl;	
        output_particles.clear();
	for (auto &c: sorted_by_pt(jet.constituents())){
	  output_particles.push_back(input_particles.at(c.user_index()));
	}
        std::cout << "output_particles.size()" << output_particles.size() << std::endl;
	output_particles.resize(NMAX);
	
	fill(vpt, output_particles, [](PFCand& p) { return p.pt; });
	fill(veta, output_particles, [](PFCand& p) { return p.eta; });
	fill(vphi, output_particles, [](PFCand& p) { return p.phi; });
	fill(ve, output_particles, [](PFCand& p) { return p.e; });
	fill(vpdgid, output_particles, [](PFCand& p) { return p.pdgid; });
	fill(vcharge, output_particles, [](PFCand& p) { return p.charge; });
	fill(vd0, output_particles, [](PFCand& p) { return p.d0; });
	fill(vdz, output_particles, [](PFCand& p) { return p.dz; });

	tout->Fill();
      }
      //if (jettype > -1) break;
      //if (dau1_parton_pt > 0.) break;
    }
    //tout->Fill();

    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
