/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class PidFilter
 *
 *  Removes all generated particles except electrons, muons, taus,
 *  and particles with status == 3.
 *
 *  \author J. Hirschauer - FNAL
 *
 */

#include "modules/PidFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

PidFilter::PidFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PidFilter::~PidFilter()
{
}

//------------------------------------------------------------------------------

void PidFilter::Init()
{
  // PT threshold
  fPTMin = GetDouble("PTMin", 0.001);

  // keep or remove pileup particles
  fRequireNotPileup = GetBool("RequireNotPileup", false);

  // target PID
  fTarget = GetInt("TargetId", 23);

  // save more than first copy of particle?
  fOnlyFirst = GetBool("OnlyFirst", true);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void PidFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void PidFilter::Process()
{
  Candidate *candidate;
  Int_t pdgCode;
  Bool_t pass;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    pdgCode = TMath::Abs(candidate->PID);

    pass = (pdgCode == fTarget);

    if(!pass || (candidate->Momentum.Pt() < fPTMin)) continue;

    // not pileup particles
    if(fRequireNotPileup && (candidate->IsPU > 0)) continue;

    fOutputArray->Add(candidate);
    if (fOnlyFirst)
      break; // Only save first one
  }
}
