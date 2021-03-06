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
//                 DetectorConstruction.cc
// --------------------------------------------------------------
//
// Code developed by:  Victor Gabriel Leandro Alves

#include "DetectorConstruction.hh"
#include <G4Material.hh>
#include <G4VisAttributes.hh>

DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::DetectorConstruction(const G4String File)
    : fGDMLFile(File) {
    // read gdml files at constructors
    readGDML();
}

DetectorConstruction::~DetectorConstruction() {
    //    if(detectorMessenger) delete detectorMessenger;
}

void DetectorConstruction::readGDML() {
    // **** LOOK HERE*** FOR READING GDML FILES
    //    fParser.Read(fGDMLFile, false);
    fParser.Read(fGDMLFile, true);
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    // Giving World Physical Volume from GDML Parser
    fWorldPhysVol = fParser.GetWorldVolume();
    fWorldPhysVol->GetLogicalVolume()->SetVisAttributes(
        G4VisAttributes::Invisible);
}

G4VPhysicalVolume *DetectorConstruction::Construct() { return fWorldPhysVol; }
