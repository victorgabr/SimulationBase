#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

RunAction::RunAction() {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run *aRun) {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void RunAction::EndOfRunAction(const G4Run *aRun) {
    G4cout << "number of events = " << aRun->GetNumberOfEvent() << G4endl;
}
