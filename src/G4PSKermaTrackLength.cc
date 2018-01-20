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
// --------------------------------------------------------------
//                 GEANT 4 - BrachySourceKerma
// --------------------------------------------------------------
//
// Code developed by:  Victor Gabriel Leandro Alves
// Copyright 2008-2017
// GEANT4 class G4PSTrackLength.hh was extended to calculate Kerma using Linear
// Kerma track-length
// estimator.
//      ref: https://www.ncbi.nlm.nih.gov/pubmed/3626996
//    *******************************
//    *                             *
//    *    G4PSKermaTrackLength.cc                *
//    *                             *
//    *******************************
#include "G4PSKermaTrackLength.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UnitsTable.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include <G4EmCalculator.hh>
#include <G4NistManager.hh>
#include <iostream>
#include <math.h>
#include <stdio.h>

G4PSKermaTrackLength::G4PSKermaTrackLength(G4String name, G4int depth)
    : G4VPrimitiveScorer(name, depth), HCID(-1), weighted(false),
      multiplyKinE(true), divideByVelocity(false) {
    ;
}

G4PSKermaTrackLength::~G4PSKermaTrackLength() { delete _matManager; }

G4double G4PSKermaTrackLength::calcMassAttenuationTotal(G4double energy,
                                                        G4Material *mat) {
    G4EmCalculator emCalculator;
    G4double tl = emCalculator.ComputeGammaAttenuationLength(energy, mat);
    // length in cm
    tl /= cm;
    G4double u_total = 1.0 / tl;
    auto density = mat->GetDensity() / (g / cm3);

    return u_total / density;
}

G4double
G4PSKermaTrackLength::calcMassEnergyAbsorptionCoefficient(G4double Z,
                                                          G4double uTotal) {
    // z 1 - 9
    return calcUenZ19(Z, uTotal);
}

G4double G4PSKermaTrackLength::calcUenZ19(G4double Z, G4double uTotal) {
    // fit coeffs ref Eur. Phys. J. D (2017) 71: 235
    // DOI: 10.1140/epjd/e2017-70679-7
    G4double a = -1.636147435e-3;
    G4double b = 1.014600795;
    G4double c = 0.8557756754;
    G4double d = 0.8245590874;
    G4double fac = (a * pow(Z, 2) + b * Z - c) / (Z - d);

    return fac * uTotal;
}

G4bool G4PSKermaTrackLength::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
    G4double trklength = aStep->GetStepLength();

    if (trklength == 0.) {
        return FALSE;
    }

    if (weighted) {
        trklength *= aStep->GetPreStepPoint()->GetWeight();
    }

    // Calculo do Kerma track length estimator

    // Obtendo os volumes dos detectores

    G4VPhysicalVolume *physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
    G4VPVParameterisation *physParam = physVol->GetParameterisation();
    G4VSolid *solid = 0;

    if (physParam) { // for parameterized volume
        G4int idx =
            ((G4TouchableHistory *)(aStep->GetPreStepPoint()->GetTouchable()))
                ->GetReplicaNumber(indexDepth);
        solid = physParam->ComputeSolid(idx, physVol);
        solid->ComputeDimensions(physParam, idx, physVol);
    } else { // for ordinary volume
        solid = physVol->GetLogicalVolume()->GetSolid();
    }

    G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
    G4double uen = getUen(energy);

    // test calculating cross section
    // managers

    //    G4Material *detMaterial = physVol->GetLogicalVolume()->GetMaterial();

    //    auto propTable = detMaterial->GetMaterialPropertiesTable();
    //    //    propTable->DumpTable();
    //    auto nistDetMaterial =
    //        _matManager->FindOrBuildMaterial(detMaterial->GetName());

    //    auto nistComp = nistDetMaterial->GetMatComponents();
    //    auto components = detMaterial->GetMatComponents();

    //    const G4String mname = detMaterial->GetName();

    //    const G4double *mf = detMaterial->GetFractionVector();
    //    size_t ne = detMaterial->GetNumberOfElements();

    //    auto mtable = detMaterial->GetMaterialTable();
    //    std::vector<G4String> namesVec;
    //    for (auto &m : *mtable) {
    //        auto n = m->GetName();
    //        namesVec.push_back(n);
    //    }

    //    std::vector<G4double> zVec;
    //    auto total{0.0};
    //    std::vector<G4double> frac;
    //    std::vector<G4double> uenVec;
    //    const G4ElementVector *elVec = detMaterial->GetElementVector();
    //    for (const G4Element *el : *elVec) {

    //        total += el->GetZ();
    //        G4int Z = static_cast<G4int>(el->GetZ());
    //        zVec.push_back(Z);
    //        auto mat = _matManager->FindOrBuildSimpleMaterial(Z);
    //        // get density in g/cm3
    //        auto uTotal_rho = calcMassAttenuationTotal(energy, mat);
    //        G4double uen_rho =
    //            calcMassEnergyAbsorptionCoefficient(el->GetZ(), uTotal_rho);
    //        uenVec.push_back(uen_rho);
    //    }

    //    for (auto z : zVec) {
    //        frac.push_back(z / total);
    //    }

    //    G4double uen_t = 0.0;
    //    size_t i = 0;
    //    for (auto u : uenVec) {
    //        uen_t += u * frac[i];
    //        i++;
    //    }

    //    G4double uen_rho = calcMassEnergyAbsorptionCoefficient(Z, uTotal);

    if (multiplyKinE) {
        trklength *= aStep->GetPreStepPoint()->GetKineticEnergy() * uen /
                     (solid->GetCubicVolume());
    }
    //----------------------------------------------------------------------------------------------//

    if (divideByVelocity) {
        trklength /= aStep->GetPreStepPoint()->GetVelocity();
    }

    const G4int index = GetIndex(aStep);
    EvtMap->add(index, trklength);
    return TRUE;
}

G4double G4PSKermaTrackLength::getUen(G4double en) {
    // Hardcoded  Uen data for water
    //  ref: https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
    // TODO FIT data
    G4double log10e = log10(en);

    // Fitted Uen/r - NIST/Water
    G4double Log10uen =
        (-1.508473990248644 +
         log10e *
             (-3.990441810810287 +
              log10e *
                  (-4.939193679834307 +
                   log10e *
                       (-3.660466557887870 +
                        log10e *
                            (-1.679146548266836 +
                             log10e * (-0.3992321053137942 +
                                       log10e * -0.02725354161828657)))))) /
        (1.0 +
         log10e *
             (2.527911872893696 +
              log10e *
                  (2.797856401815172 +
                   log10e *
                       (1.700542145298706 +
                        log10e * (0.5878067109045983 +
                                  log10e * (0.09967938698396550 +
                                            log10e * 0.007659305032379002))))));

    G4double uen = pow(10, Log10uen);
    return uen;
}

void G4PSKermaTrackLength::Initialize(G4HCofThisEvent *HCE) {
    EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if (HCID < 0) {
        HCID = GetCollectionID(0);
    }
    HCE->AddHitsCollection(HCID, (G4VHitsCollection *)EvtMap);
}

void G4PSKermaTrackLength::EndOfEvent(G4HCofThisEvent *) { ; }

void G4PSKermaTrackLength::clear() { EvtMap->clear(); }

void G4PSKermaTrackLength::DrawAll() { ; }

void G4PSKermaTrackLength::PrintAll() {
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
    G4cout << " PrimitiveScorer " << GetName() << G4endl;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl;
    std::map<G4int, G4double *>::iterator itr = EvtMap->GetMap()->begin();
    for (; itr != EvtMap->GetMap()->end(); itr++) {
        G4cout << "  copy no.: " << itr->first << "  track length: ";
        if (multiplyKinE && !divideByVelocity) {
            G4cout << *(itr->second) / (mm * MeV) << " mm*MeV";
        } else if (!multiplyKinE && divideByVelocity) {
            G4cout << *(itr->second) * second << " /second";
        } else if (multiplyKinE && divideByVelocity) {
            G4cout << *(itr->second) / MeV * second << " MeV/second";
        } else {
            G4cout << G4BestUnit(*(itr->second), "Length");
        }
        G4cout << G4endl;
    }
}
