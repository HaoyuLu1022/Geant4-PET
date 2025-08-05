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
//
/// \file B3/B3b/src/DetectorConstruction.cc
/// \brief Implementation of the B3::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

namespace B3
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
G4NistManager* man = G4NistManager::Instance();

G4bool isotopes = false;

G4Element* O = man->FindOrBuildElement("O", isotopes);
G4Element* Si = man->FindOrBuildElement("Si", isotopes);
G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

auto LSO = new G4Material("Lu2SiO5", 7.4 * g / cm3, 3);
LSO->AddElement(Lu, 2);
LSO->AddElement(Si, 1);
LSO->AddElement(O, 5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
// Siemens Biograph Vision detector Parameters (60,800 total crystals)
//
G4double cryst_dX = 3.2 * mm, cryst_dY = 3.2 * mm, cryst_dZ = 20. * mm;

// Mini-block parameters (5x5 crystals per mini-block)
G4int crystals_per_miniblock_axial = 5;   // axial direction
G4int crystals_per_miniblock_trans = 5;   // transaxial direction
G4int crystals_per_miniblock = crystals_per_miniblock_axial * crystals_per_miniblock_trans; // 25 crystals per mini-block

// Block parameters (4x2 mini-blocks per block = 20x10 crystals per block)
G4int miniblocks_per_block_axial = 2;  // 2 mini-blocks axially
G4int miniblocks_per_block_trans = 4;  // 4 mini-blocks transaxially
G4int miniblocks_per_block = miniblocks_per_block_axial * miniblocks_per_block_trans; // 8 mini-blocks per block

// Scanner geometry
G4int blocks_per_ring = 38;  // number of blocks around the ring
G4double scanner_inner_radius = 41.0 * cm;  // inner radius of detector ring
G4int nb_rings = 8;  // number of axial rings

// Calculated parameters
G4double miniblock_dX = crystals_per_miniblock_axial * cryst_dX;     // mini-block axial size
G4double miniblock_dY = crystals_per_miniblock_trans * cryst_dY;     // mini-block transaxial size
G4double block_dX = miniblocks_per_block_axial * miniblock_dX;       // block axial size (20 crystals)
G4double block_dY = miniblocks_per_block_trans * miniblock_dY;       // block transaxial size (10 crystals)

// Ring geometry calculations
G4double dPhi_block = twopi / blocks_per_ring;  // angular spacing between blocks
G4double block_depth = cryst_dZ + 5.0 * mm;     // block depth (crystal + housing)
G4double ring_inner_radius = scanner_inner_radius;
G4double ring_outer_radius = scanner_inner_radius + block_depth;

// Total detector dimensions
G4double detector_dZ = nb_rings * block_dX;
//
G4NistManager* nist = G4NistManager::Instance();
G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
G4Material* cryst_mat = nist->FindOrBuildMaterial("Lu2SiO5");

//
// World
//
G4double world_sizeXY = 2.4 * ring_outer_radius;
G4double world_sizeZ = 1.2 * detector_dZ;

auto solidWorld =
    new G4Box("World",  // its name
            0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        default_mat,  // its material
                                        "World");  // its name

auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                    G4ThreeVector(),  // at (0,0,0)
                                    logicWorld,  // its logical volume
                                    "World",  // its name
                                    nullptr,  // its mother  volume
                                    false,  // no boolean operation
                                    0,  // copy number
                                    fCheckOverlaps);  // checking overlaps

//
// Define Crystal (Level 1: Individual LSO crystal)
//
G4double gap = 0.5 * mm;  // a gap for wrapping
G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
auto solidCryst = new G4Box("Crystal", dX / 2, dY / 2, cryst_dZ / 2);

auto logicCryst = new G4LogicalVolume(solidCryst,  // its solid
                                        cryst_mat,  // its material
                                        "CrystalLV");  // its name

//
// Define Mini-block (Level 2: 5x5 crystal array)
//
auto solidMiniBlock = new G4Box("MiniBlock", miniblock_dX / 2, miniblock_dY / 2, cryst_dZ / 2);

auto logicMiniBlock = new G4LogicalVolume(solidMiniBlock,  // its solid
                                          default_mat,  // its material (air)
                                          "MiniBlockLV");  // its name

// Place crystals within a mini-block (5x5 arrangement)
for (G4int ix = 0; ix < crystals_per_miniblock_axial; ix++) {
    for (G4int iy = 0; iy < crystals_per_miniblock_trans; iy++) {
        G4double x_pos = (ix - (crystals_per_miniblock_axial - 1) / 2.0) * cryst_dX;
        G4double y_pos = (iy - (crystals_per_miniblock_trans - 1) / 2.0) * cryst_dY;
        G4int copy_number = ix * crystals_per_miniblock_trans + iy;
        
        new G4PVPlacement(nullptr,  // no rotation
                        G4ThreeVector(x_pos, y_pos, 0),  // position
                        logicCryst,  // its logical volume
                        "Crystal",  // its name
                        logicMiniBlock,  // its mother volume
                        false,  // no boolean operation
                        copy_number,  // copy number
                        fCheckOverlaps);  // checking overlaps
    }
}

//
// Define Block (Level 3: 4x2 mini-block array = 20x10 crystals)
//
auto solidBlock = new G4Box("Block", block_dX / 2, block_dY / 2, cryst_dZ / 2);

auto logicBlock = new G4LogicalVolume(solidBlock,  // its solid
                                      default_mat,  // its material (air)
                                      "BlockLV");  // its name

// Place mini-blocks within a block (4x2 arrangement)
for (G4int ix = 0; ix < miniblocks_per_block_axial; ix++) {
    for (G4int iy = 0; iy < miniblocks_per_block_trans; iy++) {
        G4double x_pos = (ix - (miniblocks_per_block_axial - 1) / 2.0) * miniblock_dX;
        G4double y_pos = (iy - (miniblocks_per_block_trans - 1) / 2.0) * miniblock_dY;
        G4int copy_number = ix * miniblocks_per_block_trans + iy;
        
        new G4PVPlacement(nullptr,  // no rotation
                        G4ThreeVector(x_pos, y_pos, 0),  // position
                        logicMiniBlock,  // its logical volume
                        "MiniBlock",  // its name
                        logicBlock,  // its mother volume
                        false,  // no boolean operation
                        copy_number,  // copy number
                        fCheckOverlaps);  // checking overlaps
    }
}

//
// Define Ring (Level 4: Ring of blocks)
//
auto solidRing = new G4Tubs("Ring", ring_inner_radius, ring_outer_radius, block_dX / 2, 0., twopi);

auto logicRing = new G4LogicalVolume(solidRing,  // its solid
                                    default_mat,  // its material
                                    "RingLV");  // its name

// Place blocks within a ring
for (G4int iblock = 0; iblock < blocks_per_ring; iblock++) {
    G4double phi = iblock * dPhi_block;
    G4RotationMatrix rotm = G4RotationMatrix();
    rotm.rotateY(90 * deg);
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0.);
    G4ThreeVector position = (ring_inner_radius + block_depth / 2) * uz;
    G4Transform3D transform = G4Transform3D(rotm, position);

    new G4PVPlacement(transform,  // rotation,position
                    logicBlock,  // its logical volume
                    "Block",  // its name
                    logicRing,  // its mother volume
                    false,  // no boolean operation
                    iblock,  // copy number
                    fCheckOverlaps);  // checking overlaps
}

//
// Full detector (Level 5: Multiple rings)
//
auto solidDetector = new G4Tubs("Detector", ring_inner_radius, ring_outer_radius, 0.5 * detector_dZ, 0., twopi);

auto logicDetector = new G4LogicalVolume(solidDetector,  // its solid
                                        default_mat,  // its material
                                        "DetectorLV");  // its name

//
// place rings within detector
//
G4double ring_spacing = block_dX;
G4double first_ring_z = -0.5 * (detector_dZ - ring_spacing);
for (G4int iring = 0; iring < nb_rings; iring++) {
    G4double ring_z = first_ring_z + iring * ring_spacing;
    new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(0, 0, ring_z),  // position
                    logicRing,  // its logical volume
                    "Ring",  // its name
                    logicDetector,  // its mother volume
                    false,  // no boolean operation
                    iring,  // copy number
                    fCheckOverlaps);  // checking overlaps
}

//
// place detector in world
//
new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicDetector,  // its logical volume
                    "Detector",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    fCheckOverlaps);  // checking overlaps

//
// patient
//
G4double patient_radius = 8 * cm;
G4double patient_dZ = 10 * cm;
G4Material* patient_mat = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");

auto solidPatient = new G4Tubs("Patient", 0., patient_radius, 0.5 * patient_dZ, 0., twopi);

auto logicPatient = new G4LogicalVolume(solidPatient,  // its solid
                                        patient_mat,  // its material
                                        "PatientLV");  // its name

//
// place patient in world
//
new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicPatient,  // its logical volume
                    "Patient",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    fCheckOverlaps);  // checking overlaps

// Visualization attributes
//
logicMiniBlock->SetVisAttributes(G4VisAttributes::GetInvisible());
logicBlock->SetVisAttributes(G4VisAttributes::GetInvisible());
logicRing->SetVisAttributes(G4VisAttributes::GetInvisible());
logicDetector->SetVisAttributes(G4VisAttributes::GetInvisible());

// Set crystal color to blue for visibility
auto crystalVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.5)); // Semi-transparent blue
logicCryst->SetVisAttributes(crystalVisAtt);

// Print materials
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

// always return the physical World
//
return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

// declare crystal as a MultiFunctionalDetector scorer
//
auto cryst = new G4MultiFunctionalDetector("crystal");
G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
cryst->RegisterPrimitive(primitiv1);
SetSensitiveDetector("CrystalLV", cryst);

// declare patient as a MultiFunctionalDetector scorer
//
auto patient = new G4MultiFunctionalDetector("patient");
G4SDManager::GetSDMpointer()->AddNewDetector(patient);
G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
patient->RegisterPrimitive(primitiv2);
SetSensitiveDetector("PatientLV", patient);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B3
