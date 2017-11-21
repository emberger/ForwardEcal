#include "B4cReadoutGeometry.hh"
#include "B4cCalorimeterSD.hh"


#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include <iostream>

#include "B4cDetParams.hh"

MyRO::MyRO() : G4VReadOutGeometry(){
}

MyRO::MyRO(G4String aString) : G4VReadOutGeometry(aString){



}

MyRO::~MyRO(){
}





G4VPhysicalVolume* MyRO::Build(){

								// Geometry parameters also for Detector Construction
								GetInst().SetfNofLayers(100);		// in #
								GetInst().SetcalorSizeXY(2000); // in mm
								GetInst().SettileLenX(20); // in mm
								GetInst().SettileLenY(20); // in mm

								GetInst().SetabsoThickness(1); // in mm
								GetInst().SetgapThickness(5); // in mm

								GetInst().SetAbsFirst(true);    // true if structure should be |Abso-Gap|Abso-Gap|... false if structure should be |Gap-Abso|Gap-Abso|...


								GetInst().SetWorldMult(10.);//in mm
								GetInst().InitDet(); // dont forget!!!




//	   auto ROworldSizeXY = 1.2 * GetInst().GetcalorSizeXY();
//	    auto ROworldSizeZ  = 1.2 * GetInst().GetcalorThickness();

								auto fCheckOverlaps=true;

								// A dummy material is used to fill the volumes of the readout geometry.
								// ( It will be allowed to set a NULL pointer in volumes of such virtual
								// division in future, since this material is irrelevant for tracking.)
								G4Material * dummyMat  = new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);



								//Build Readout World

								auto ROWorldS = new G4Box("ROWorld",
																																		GetInst().GetWorldSizeXY()/2,
																																		GetInst().GetWorldSizeXY()/2,
																																		GetInst().GetWorldSizeZ()/2);

								auto ROWorldLog = new G4LogicalVolume(ROWorldS,
																																														dummyMat,
																																														"ROWorldLogical",
																																														0,
																																														0,
																																														0);


								G4PVPlacement *ROWorld = new G4PVPlacement(0,
																																																			G4ThreeVector(),
																																																			"ROWorldPhysical",
																																																			ROWorldLog,
																																																			0,
																																																			false,
																																																			0);


								//Place Calorimeter

								auto calorimeterS
																= new G4Box("Calorimeter", // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetcalorThickness()/2); // its size

								auto calorLV
																= new G4LogicalVolume(
																calorimeterS, // its solid
																dummyMat, // its material
																"Calorimeter"); // its name

								G4ThreeVector seg1(0*mm,0*mm,(GetInst().GetcalorThickness()/2)*mm);
								new G4PVPlacement(
																0,                             // no rotation
																seg1,               // at (0,0,0)
																calorLV,                       // its logical volume
																"Calorimeter",                 // its name
																ROWorldLog,                       // its mother  volume
																false,                         // no boolean operation
																1,                             // copy number
																fCheckOverlaps);            // checking overlaps





								// Build Layers

								auto ROlayerS
																= new G4Box("Layer",  // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetlayerThickness()/2); //its size

								auto ROlayerLV
																= new G4LogicalVolume(
																ROlayerS,               // its solid
																dummyMat,      // its material
																"Layer");             // its name

								new G4PVReplica(
																"Layer",              // its name
																ROlayerLV,              // its logical volume
																calorLV,              // its mother
																kZAxis,               // axis of replication
																GetInst().GetfNofLayers(),            // number of replica
																GetInst().GetlayerThickness());      // witdth of replica

								//-------------------------------
								//build calorimeter readout cells
								//-------------------------------

								//build Gap
								auto GapS
																=new G4Box("Gap", GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetgapThickness()/2);

								auto GapLV
																=new G4LogicalVolume(GapS,
																																					dummyMat,
																																					"Gap");
								if(GetInst().GetAbsFirst()) {

																new G4PVPlacement(0,
																																		G4ThreeVector(0., 0., GetInst().GetabsoThickness()/2 ),
																																		GapLV,
																																		"Gap",
																																		ROlayerLV,
																																		false,
																																		0,
																																		true);


								}
								else{
																new G4PVPlacement(0,
																																		G4ThreeVector(0., 0., -GetInst().GetabsoThickness()/2 ),
																																		GapLV,
																																		"Gap",
																																		ROlayerLV,
																																		false,
																																		0,
																																		true);

								}
								//build strips parallel to X axis

								auto StripS
																=new G4Box("Strip", GetInst().GetcalorSizeXY()/2,GetInst().GettileLenY()/2,GetInst().GetgapThickness()/2);

								auto StripLV
																=new G4LogicalVolume(StripS,
																																					dummyMat,
																																					"Strip");

								new G4PVReplica("Strip",
																								StripLV, //logical volume
																								GapLV, //mother volume
																								kYAxis,
																								GetInst().GetnofTilesY(),
																								GetInst().GettileLenY());


								// build cells in Xaxis strip

								auto CellS
																=new G4Box("Cell", GetInst().GettileLenX()/2, GetInst().GettileLenY()/2, GetInst().GetgapThickness()/2);

								auto CellLV
																=new G4LogicalVolume(CellS,
																																					dummyMat,
																																					"Cell");

								new G4PVReplica("Cell",
																								CellLV, //logical volume
																								StripLV, //mother volume
																								kXAxis,
																								GetInst().GetnofTilesX(),
																								GetInst().GettileLenX());


								auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
								simpleBoxVisAtt->SetVisibility(true);
								CellLV->SetVisAttributes(simpleBoxVisAtt);



								auto dummy
																= new B4cCalorimeterSD("dummy", "dummyCollection",0,0,0);

								CellLV->SetSensitiveDetector(dummy);



								return ROWorld;




}
