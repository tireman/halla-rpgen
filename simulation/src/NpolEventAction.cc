//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

//%% NpolEventAction.cc %% 

// This is were the event level actions are placed
// Created: Daniel Wilbern November 2014
// Modified: William Tireman December 2014

#include <cstring>
#include <string>
#include "G4Event.hh"
#include "G4RunManager.hh"

#include "NpolEventAction.hh"
#include "NpolAnalysisManager.hh"

NpolEventAction::NpolEventAction()
{
  // set printing per each 10 events default (can be changed in macro)
  G4RunManager::GetRunManager()->SetPrintProgress(10);   
}

NpolEventAction::~NpolEventAction()
{}

void NpolEventAction::BeginOfEventAction(const G4Event* evt) {
  NpolAnalysisManager::GetInstance()->BeginEvent(evt->GetEventID());
}

void NpolEventAction::EndOfEventAction(const G4Event* evt) {
  NpolAnalysisManager::GetInstance()->EndEvent(evt->GetEventID());
}

