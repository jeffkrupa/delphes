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

using namespace std;


//---------------------------------------------------------------------------


static int NMAX = 9000;

struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float puppi = 1;
  float pdgid = 0;
  float hardfrac = 1;  
  float cluster_idx = -1;
  float cluster_hardch_pt = 0;
  float cluster_puch_pt = 0;
  float cluster_r = 0;
  float vtxid = -1;
  float npv = 0;
};


class Cluster : public vector<PFCand*> {
public:
  void finalize () 
  {
    for (auto* p : *this) {
      _sum_pt += p->pt;
      if (p->vtxid==0){
	_hardch_pt += p->pt;
	//std::cout << _hardch_pt << std::endl;
      }
      else if (p->vtxid==1){
	_puch_pt += p->pt;
      }
      _eta += p->pt * p->eta;
      _phi += p->pt * p->phi;
    }
    _eta /= _sum_pt;
    _phi /= _sum_pt;

    //if (_hardch_pt>90)
    //std::cout << _hardch_pt << std::endl;

    float largest_dr = -1.;
    for (auto* p : *this) {
      auto dr = pow(_eta - p->eta, 2) + pow(_phi - p->phi, 2);
      if (dr > largest_dr)
	largest_dr = dr;
    }
    _r = largest_dr;
  }

  float eta() { return _eta; }
  float phi() { return _phi; }
  float sum_pt() { return _sum_pt; }
  float hardch_pt() { return _hardch_pt; }
  float puch_pt() { return _puch_pt; }
  float r() { return _r; }

private:
  float _eta=0, _phi=0, _sum_pt=0, _hardch_pt=0, _puch_pt=0, _r=0;

};


template<int K>
class KMeans { 
public:
    KMeans(vector<PFCand*> particles, int max_iter=20) 
    {
        // randomly initialize centroids
        array<int, K> i_centroids;
        for (int i=0; i!=K; ++i) {
            while (true) {
                int i_p = rand() % particles.size();
                bool found = false;
                for (int j=0; j!=i; ++j) {
                    found = (i_centroids[j] == i_p);
                    if (found)
                        break;
                }
                if (!found) {
                    i_centroids[i] = i_p;
                    centroids[i][0] = particles[i_p]->eta;
                    centroids[i][1] = particles[i_p]->phi;
                    break;
                }
            }
        } 

        for (int i_iter=0; i_iter!=max_iter; ++i_iter) {
            assign_particles(particles);
            update_centroids();
        }
    }

    ~KMeans() { }

    const array<Cluster, K> get_clusters() { return clusters; }

private:
    array<array<float, 2>, K> centroids;
    array<Cluster, K> clusters;

    void assign_particles(vector<PFCand*> &particles) 
    {
        for (int i=0; i!=K; ++i) {
            clusters[i].clear();
        }

        for (auto& p : particles) {
            float closest = 99999;
            int i_closest = -1;
            float eta = p->eta; float phi = p->phi;

            for (int i=0; i!=K; ++i) {
                auto dr = pow(eta - centroids[i][0], 2) + pow(phi - centroids[i][1], 2);
                if (dr < closest) {
                    closest = dr;
                    i_closest = i;
                }
            }
            clusters[i_closest].push_back(p);
        }
    }

    void update_centroids() 
    {
        for (int i=0; i!=K; ++i) {
            float eta_sum=0, phi_sum=0;
            auto &cluster = clusters[i];
            for (auto& p : cluster) {
                eta_sum += p->eta;
                phi_sum += p->phi;
            }
            centroids[i][0] = eta_sum / cluster.size();
            centroids[i][1] = phi_sum / cluster.size();
        }
    }
}; 




template <int K, int N>
class HierarchicalOrdering {
public:
    HierarchicalOrdering() { }

    vector<Cluster> 
    fit(vector<PFCand> &particles)
    {
        vector<PFCand*> p_particles;
        for (auto& p : particles)
            p_particles.push_back(&p);

        auto clusters = _recursive_fit(p_particles);
        for (auto& c : clusters)
          c.finalize();
        return clusters;
    }
    
private:
    vector<Cluster> 
    _recursive_fit(const vector<PFCand*> &particles)
    {
        vector<Cluster> clusters;

        auto kmeans = KMeans<K>(particles);
        for (int i_k=0; i_k!=K; ++i_k) {
            auto cluster = kmeans.get_clusters()[i_k];
            if (cluster.size() > N) {
                auto split_clusters = _recursive_fit(cluster);
                for (auto& c : split_clusters) {
                    clusters.push_back(c);
                }
            } else {
                clusters.push_back(cluster);
            } 
        }
        return clusters;
    }
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
  TBranch* pfbranch = (TBranch*)itree->GetBranch("ParticleFlowCandidate");
  TBranch* genjetbranch = (TBranch*)itree->GetBranch("GenJet");
  std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;

  vector<PFCand> output_particles;
  output_particles.reserve(NMAX);

  vector<float> vpt, veta, vphi, ve, vpuppi, vpdgid, vhardfrac, vcluster_idx, vvtxid, vcluster_r, vcluster_hardch_pt, vcluster_puch_pt, vnpv;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX); vpuppi.reserve(NMAX); vpdgid.reserve(NMAX); 
  vhardfrac.reserve(NMAX); vcluster_idx.reserve(NMAX); vvtxid.reserve(NMAX);
  vcluster_r.reserve(NMAX); vcluster_hardch_pt.reserve(NMAX); vcluster_puch_pt.reserve(NMAX);
  vnpv.reserve(NMAX);

  float genmet=-99., genmetphi=-99.;
  float genjet1pt=-99., genjet1eta=-99., genjet1phi=-99., genjet1e=-99.;
  float genjet2pt=-99., genjet2eta=-99., genjet2phi=-99., genjet2e=-99.;
  
  tout->Branch("pt", &vpt);
  tout->Branch("eta", &veta);
  tout->Branch("phi", &vphi);
  tout->Branch("e", &ve);
  tout->Branch("puppi", &vpuppi);
  tout->Branch("pdgid", &vpdgid);
  tout->Branch("hardfrac", &vhardfrac);
  tout->Branch("cluster_idx", &vcluster_idx);
  tout->Branch("cluster_r", &vcluster_r);
  tout->Branch("cluster_hardch_pt", &vcluster_hardch_pt);
  tout->Branch("cluster_puch_pt", &vcluster_puch_pt);
  tout->Branch("vtxid", &vvtxid);
  tout->Branch("npv", &vnpv);
  TBranch* b_genmet = tout->Branch("genmet",&genmet, "genmet/F");
  TBranch* b_genmetphi = tout->Branch("genmetphi",&genmetphi, "genmetphi/F");
  TBranch* b_genjet1pt = tout->Branch("genjet1pt",&genjet1pt, "genjet1pt/F");
  TBranch* b_genjet1eta = tout->Branch("genjet1eta",&genjet1eta, "genjet1eta/F");
  TBranch* b_genjet1phi = tout->Branch("genjet1phi",&genjet1phi, "genjet1phi/F");
  TBranch* b_genjet1e = tout->Branch("genjet1e",&genjet1e, "genjet1e/F");
  TBranch* b_genjet2pt = tout->Branch("genjet2pt",&genjet2pt, "genjet2pt/F");
  TBranch* b_genjet2eta = tout->Branch("genjet2eta",&genjet2eta, "genjet2eta/F");
  TBranch* b_genjet2phi = tout->Branch("genjet2phi",&genjet2phi, "genjet2phi/F");
  TBranch* b_genjet2e = tout->Branch("genjet2e",&genjet2e, "genjet2e/F");

  auto ho = HierarchicalOrdering<4, 10>();
  //auto ho = HierarchicalOrdering<4, 30>();

  ExRootProgressBar progressBar(nevt);
  
  auto comp_pt = [](auto &a, auto &b) { return a.sum_pt() > b.sum_pt(); };
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);
    
    float npv = itree->GetLeaf("Vertex_size")->GetValue(0);
    genmet = itree->GetLeaf("GenMissingET.MET")->GetValue(0);
    genmetphi = itree->GetLeaf("GenMissingET.Phi")->GetValue(0);

    unsigned int ngenjets = genjetbranch->GetEntries();
    ngenjets = itree->GetLeaf("GenJet_size")->GetValue(0);
    for (unsigned int j=0; j<ngenjets; j++){
      if (j>1)
	break;
      TLorentzVector tmpjet;
      tmpjet.SetPtEtaPhiM(itree->GetLeaf("GenJet.PT")->GetValue(j),itree->GetLeaf("GenJet.Eta")->GetValue(j),itree->GetLeaf("GenJet.Phi")->GetValue(j),itree->GetLeaf("GenJet.Mass")->GetValue(j));
      if (j==0){
	genjet1pt = tmpjet.Pt();
	genjet1eta = tmpjet.Eta();
	genjet1phi = tmpjet.Phi();
	genjet1e = tmpjet.E();
      }
      if (j==1){
	genjet2pt = tmpjet.Pt();
	genjet2eta = tmpjet.Eta();
	genjet2phi = tmpjet.Phi();
	genjet2e = tmpjet.E();
      }      
    }

    input_particles.clear();
    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);
    for (unsigned int j=0; j<npfs; j++){
      PFCand tmppf;
      tmppf.npv = npv;
      tmppf.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      tmppf.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmppf.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmppf.e = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.puppi = itree->GetLeaf("ParticleFlowCandidate.PuppiW")->GetValue(j);
      tmppf.hardfrac = itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j);
      tmppf.pdgid = itree->GetLeaf("ParticleFlowCandidate.PID")->GetValue(j);
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

    // sorting input particles by pT
    sort(input_particles.begin(), input_particles.end(), comp_p4);    

    // get clusters of 10 particles
    auto clusters = ho.fit(input_particles);

    // sort clusters by sum pT. not very efficient since we recompute
    // sum_pt for each comparison but whatever 
    sort(clusters.begin(), clusters.end(), comp_pt);

    // now sort clusters by proximity to last cluster, starting with
    // the hardest cluster
    vector<Cluster> sorted_clusters = {clusters[0]};
    clusters.erase(clusters.begin());
    while (clusters.size() > 0) {
      unsigned i_cluster = 0;
      float best_dr2 = 999999;
      auto &last_cluster = sorted_clusters.back();
      for (unsigned j=0; j!=clusters.size(); ++j) {
        float dr2 = pow(last_cluster.eta() - clusters[j].eta(), 2)
                    + pow(last_cluster.phi() - clusters[j].phi(), 2);
        if (dr2 < best_dr2) {
          i_cluster = j;
          best_dr2 = dr2;
        }
      }
      sorted_clusters.push_back(clusters[i_cluster]);
      clusters.erase(clusters.begin()+i_cluster);
    }

    output_particles.clear();
    int cluster_idx = 0;
    for (auto& cluster : sorted_clusters) {
      for (auto* p : cluster) {
        p->cluster_idx = cluster_idx;
	p->cluster_hardch_pt = cluster.hardch_pt();
	p->cluster_puch_pt = cluster.puch_pt();
	p->cluster_r = cluster.r();

	//if (p->cluster_hardch_pt>90)
	//std::cout << "WTF: " << p->cluster_hardch_pt << std::endl;
	//std::cout << "Hard: " << p->cluster_hardch_pt << std::endl;
	//std::cout << "PU: " << p->cluster_puch_pt << std::endl;
	
        output_particles.push_back(*p); 
      }
      ++cluster_idx;
    }
    // if there are fewer than NMAX, it'll get padded out with default values
    output_particles.resize(NMAX);

    fill(vpt, output_particles, [](PFCand& p) { return p.pt; }); 
    fill(veta, output_particles, [](PFCand& p) { return p.eta; }); 
    fill(vphi, output_particles, [](PFCand& p) { return p.phi; }); 
    fill(ve, output_particles, [](PFCand& p) { return p.e; }); 
    fill(vpuppi, output_particles, [](PFCand& p) { return p.puppi; }); 
    fill(vpdgid, output_particles, [](PFCand& p) { return p.pdgid; }); 
    fill(vhardfrac, output_particles, [](PFCand& p) { return p.hardfrac; }); 
    fill(vcluster_idx, output_particles, [](PFCand& p) { return p.cluster_idx; }); 
    fill(vcluster_r, output_particles, [](PFCand& p) { return p.cluster_r; }); 
    fill(vcluster_hardch_pt, output_particles, [](PFCand& p) { return p.cluster_hardch_pt; }); 
    fill(vcluster_puch_pt, output_particles, [](PFCand& p) { return p.cluster_puch_pt; }); 
    fill(vvtxid, output_particles, [](PFCand& p) { return p.vtxid; }); 
    fill(vnpv, output_particles, [](PFCand& p) { return p.npv; }); 
    
    tout->Fill();

    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
