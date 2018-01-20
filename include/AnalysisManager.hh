//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/*
Based on thee class BrachyAnalysisManager creates and manages histograms and
ntuples by Susanna Guatelli
*/

#ifndef AnalysisManager_HH
#define AnalysisManager_HH

#include "globals.hh"

#ifdef ANALYSIS_USE
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#endif

class AnalysisManager {
  private:
    AnalysisManager();

  public:
    ~AnalysisManager();
    static AnalysisManager *GetInstance();

    void book();
// Create the output ROOT file
// Create the ntuple and histograms

#ifdef ANALYSIS_USE
    // Method to fill the 2D histogram with the energy deposition, integrated
    // over a run, in each voxel
    // of the scoring mesh. The scoring mesh is in the plane containing the
    // source.
    void FillH2WithEnergyDeposition(G4double xx, G4double yy,
                                    G4double energyDep);

    // Energy spectrum of primary particles
    void FillPrimaryParticleHistogram(G4double);

#endif

    void save();
    // This method if called at the end of the run to store the
    // results in the ROOT file

  private:
    static AnalysisManager *instance;

#ifdef ANALYSIS_USE
    TFile *theTFile;
    TH1F *histo;
    TH2F *histo2;
#endif
};
#endif
