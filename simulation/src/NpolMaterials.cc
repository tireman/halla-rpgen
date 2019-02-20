//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% NpolMaterials.cc  %%

// Npol materials are created in here
// Created: Daniel Wilbern - November 2014
// Modified: William Tireman - December 2014-January 2015

#include <map>
#include <string>

#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "NpolMaterials.hh"

NpolMaterials *NpolMaterials::pInstance = NULL;

NpolMaterials *NpolMaterials::GetInstance() {
  if(pInstance == NULL)
	  pInstance = new NpolMaterials();
  return pInstance;
}

NpolMaterials::NpolMaterials() {
  G4cout << "Constructing NpolMaterials singleton" << G4endl;

  nistMan = G4NistManager::Instance();
  
  CreateMaterials();
  G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

NpolMaterials::~NpolMaterials() {
  std::map<std::string,G4Material *>::iterator it;
  for(it = materials.begin(); it != materials.end(); it++) 
    delete it->second;

  materials.clear();
}

G4Material *NpolMaterials::GetMaterial(const G4String material) {
  G4Material *ret = NULL;
  G4Material *mat = NULL;

  std::map<std::string,G4Material *>::iterator it = materials.find(material);
  if(it != materials.end())
  	ret = it->second;
  else if((mat = nistMan->FindOrBuildMaterial(material)) != NULL)
	  ret = mat;
  else if((mat = G4Material::GetMaterial(material)) != NULL)
	  ret = mat;
  else
    G4cout << "Material " << material << " not found." << G4endl;

  return ret;  
}

void NpolMaterials::CreateMaterials() {
  
  materials["Vacuum"] = CreateVacuum();
  materials["HardVacuum"] = CreateHardVacuum();
  materials["Air"] = CreateAir();
  materials["Scint"] = CreateScint();
  materials["Al"] = CreateAl();
  materials["LH2"] = CreateLH2();
  materials["LD2"] = CreateLD2();
  materials["Concrete"] = CreateConcrete();
  materials["Fe"] = CreateFe();
  materials["SSteel"] = CreateSSteel();
  materials["Cu"] = CreateCu();
  materials["Pb"] = CreatePb();
  materials["Kapton"] = CreateKapton();
  materials["NpolGas"] = CreateNpolGas();
  materials["Mylar"] = CreateMylar();
  materials["C2Gas"] = CreateC2Gas();
}

G4Material *NpolMaterials::CreateKapton(){

  G4String name = "";
  G4int ncomponents, natoms;
  G4Element* elH = nistMan->FindOrBuildElement("H");
  G4Element* elC = nistMan->FindOrBuildElement("C");
  G4Element* elO = nistMan->FindOrBuildElement("O");
  G4Element* elN = nistMan->FindOrBuildElement("N");
  G4double density = 1.42*g/cm3;
  G4Material *Kapton = new G4Material(name="Kapton", density, ncomponents=4);
  Kapton->AddElement(elH, natoms=10);
  Kapton->AddElement(elC, natoms=22);
  Kapton->AddElement(elO, natoms=5);
  Kapton->AddElement(elN, natoms=2);

  return Kapton;
}

G4Material *NpolMaterials::CreateNpolGas(){

  G4String name = "";
  G4double fractionmass;
  G4int ncomponents;
  G4Material* Argon = nistMan->FindOrBuildMaterial("G4_Ar");
  G4Material* CO2 = nistMan->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4double density = 1.810*mg/cm3;
  G4Material *NpolGas = new G4Material(name="NpolGas", density, ncomponents=2);
  NpolGas->AddMaterial(Argon, fractionmass=0.5);
  NpolGas->AddMaterial(CO2, fractionmass=0.5);

  return NpolGas;
}

G4Material *NpolMaterials::CreateMylar(){

  G4String name = "";
  G4int ncomponents, natoms;
  G4Element* elH = nistMan->FindOrBuildElement("H");
  G4Element* elC = nistMan->FindOrBuildElement("C");
  G4double density = 1.39*g/cm3;
  G4Material *Mylar = new G4Material(name="Mylar", density, ncomponents=2);
  Mylar->AddElement(elH, natoms=8);
  Mylar->AddElement(elC, natoms=10);

  return Mylar;
}

G4Material *NpolMaterials::CreateC2Gas(){

  G4String name = "";
  G4double fractionmass;
  G4int ncomponents;
  G4Material* Argon = nistMan->FindOrBuildMaterial("G4_Ar");
  G4Material* Ethane = nistMan->FindOrBuildMaterial("G4_ETHANE");
  G4double density = 1.497*mg/cm3;
  G4Material *C2Gas = new G4Material(name="C2Gas", density, ncomponents=2);
  C2Gas->AddMaterial(Argon, fractionmass=0.5);
  C2Gas->AddMaterial(Ethane, fractionmass=0.5);
  
  return C2Gas;
}

G4Material *NpolMaterials::CreateVacuum() {
  // Define Vacuum: New version; use Air and just make it very low density
  G4double fractionmass, density;
  G4int ncomponents;
  G4Material *Vacuum = 
	new G4Material("Vacuum", density= 1.0e-3*g/cm3, ncomponents=1, kStateGas, 0.1*kelvin, 1.e-19*atmosphere);
  Vacuum->AddMaterial(nistMan->FindOrBuildMaterial("G4_AIR"), fractionmass=1.);
  
  return Vacuum;
}

G4Material *NpolMaterials::CreateHardVacuum() {
  // Define HardVacuum: New version; use Air and just make it very low density
  G4double fractionmass, density;
  G4int ncomponents;
  G4Material *HardVacuum = new G4Material("HardVacuum", density= 1.0e-15*g/cm3,
	   ncomponents=1, kStateGas, 0.1*kelvin, 1.e-19*atmosphere);
  HardVacuum->AddMaterial(nistMan->FindOrBuildMaterial("G4_AIR"), fractionmass=1.);
  
  return HardVacuum;
}

G4Material *NpolMaterials::CreateAir() {
  return nistMan->FindOrBuildMaterial("G4_AIR");
}

G4Material *NpolMaterials::CreateAl(){
  return nistMan->FindOrBuildMaterial("G4_Al");
}

G4Material *NpolMaterials::CreateScint() {
  
  G4Element* H = nistMan->FindOrBuildElement("H");
  //G4Element* C = nistMan->FindOrBuildElement("C");
  G4Material* C = nistMan->FindOrBuildMaterial("G4_C");
  G4Material* scint = new G4Material("Scint", 1.02*g/cm3, 2);
  scint->AddElement(H, 8.451*perCent);
  //scint->AddElement(C, 91.549*perCent);
  scint->AddMaterial(C, 91.549*perCent);
  return scint;
}

G4Material *NpolMaterials::CreateLH2() {
  
  // liquid hydrogen
  // Hydrogen element  
  G4Element* H = nistMan->FindOrBuildElement("H");
  
  // Liquid Hydrogen
  G4Material *LH2 = new G4Material("LH2", 0.07085*g/cm3, 1, kStateLiquid, 15.0*kelvin);
  LH2->AddElement(H, 2);
  
  return LH2;
}

G4Material *NpolMaterials::CreateLD2() {
  // liquid deuterium
  G4int Z, N;
  G4double a;
  // Deuteron isotope
  G4Isotope* deuteron = new G4Isotope("deuteron", Z=1, N=2, a=2.0141018*g/mole);
  
  // Deuterium element
  G4Element* deuterium = new G4Element("deuterium", "deuterium", 1);
  deuterium->AddIsotope(deuteron, 1);
  
  // Liquid Deuterium
  G4Material *LD2 = new G4Material("LD2", 0.169*g/cm3, 1, kStateLiquid, 22.0*kelvin);
  LD2->AddElement(deuterium, 2);
  
  return LD2;
}

G4Material *NpolMaterials::CreateConcrete() {
  
  return nistMan->FindOrBuildMaterial("G4_CONCRETE");
}

G4Material *NpolMaterials::CreateFe() {
  
  return nistMan->FindOrBuildMaterial("G4_Fe");
}

G4Material *NpolMaterials::CreateSSteel() {
  
  return nistMan->FindOrBuildMaterial("G4_STAINLESS-STEEL");
}

G4Material *NpolMaterials::CreateCu() {
  
  return nistMan->FindOrBuildMaterial("G4_Cu");
}

G4Material *NpolMaterials::CreatePb() {
  
  return nistMan->FindOrBuildMaterial("G4_Pb");
}
