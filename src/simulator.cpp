/* 
 * fossilgen version 0.1b source code (https://github.com/trayc7/FossilGen)
 * Copyright 20013-2013
 * Tracy Heath(1,2,3) 
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 *
 * fossilgen is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck and Fredrik Ronquist
 *
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include "MbRandom.h"
#include "simulator.h"
#include "tree.h"

using namespace std;

// =========================================================================
//	Simulator::Simulator
// =========================================================================
Simulator::Simulator(MbRandom *rndp, int mt, int nt, double sc, double br, double dr, 
					 int ba, int da, double bp, double dp, bool inRt, bool allF){
	
	rando = rndp;
	modeltype = mt;
	extantStop = nt;
	treescale = sc;
	stBirthRt = br;
	stDeathRt = dr;
	birthSampAlpha = ba;
	deathSampAlpha = da;
	bPriorAlpha = bp;
	dPriorAlpha = dp;
	currentSimTime = 0.0;
	simt = NULL;
	extantTree = NULL;
	fossilTree = NULL;
	includeRootFossils = inRt;
	allFossilsSamp = allF;
	//initializeSimulator();
}

// =========================================================================
//	Simulator::Simulator
// =========================================================================
Simulator::~Simulator(){
	
	if(simt != NULL)
		delete simt;
	if(extantTree != NULL)
		delete extantTree;
	if(fossilTree != NULL)
		delete fossilTree;
}



// =========================================================================
//	Simulator::simulateTree
// =========================================================================
bool Simulator::simulateTree(){
	
	//initializeSimulator();
	bool good = false;
	while(!good){
		good = gsaBDSim();
	}
	//simt->setWholeTreeExtantRoot();

	//processSimulation();
	return good;	
}


// =========================================================================
//	Simulator::initializeSimulator
// =========================================================================
void Simulator::initializeSimulator(){
	
	Tree *phylo;
	if(simt != NULL)
		delete simt;
	if(extantTree != NULL)
		delete extantTree;

	phylo = new Tree(rando, stBirthRt, stDeathRt, extantStop, currentSimTime,
					 birthSampAlpha, deathSampAlpha, bPriorAlpha, dPriorAlpha,
					 includeRootFossils, allFossilsSamp);
	simt = phylo;
	
	

}


// =========================================================================
//	Simulator::gsaBDSim
// =========================================================================
bool Simulator::gsaBDSim(){
	
	bool treeComplete = false;
	initializeSimulator();
	simt->setNSampT(gsaExtStop);
	double dt;
	while(czechStop()){							
		dt = simt->getTimeToNextERMEvent();
		currentSimTime += dt;
		simt->ermEvent(currentSimTime);
		if(simt->getNumExtantLineages() < 1){
			treeComplete = false;
			return treeComplete;
		}
		else if(simt->getNumExtantLineages() == gsaExtStop){
			double timeIntv = simt->getTimeToNextERMEvent();
			double extT = rando->uniformRv(0, timeIntv) + currentSimTime;
			simt->setExtantPresentTime(extT);
			processGSASimulation();
		}
	}
	
	int gsaRTID = rando->discreteUniformRv(0, (int)gsaTreess.size() - 1);
	simt = gsaTreess[gsaRTID];
	simt->popExtantVec();
	extantStop = gsaExtStop;
	processSimulation(true);
	treeComplete = true;
	simt->setWholeTreeExtantRoot();
	return treeComplete;	
}


// =========================================================================
//	Simulator::czechStop
//  The Czech Stop is a Czech bakery on I-35 in West, Texas 
//  ("West" is the name of the town that is actually in North-Central Texas)
// =========================================================================
bool Simulator::czechStop(){
	
	bool keepSimulating = true;
	
	if(extantStop > 0 && simt->getNumExtantLineages() >= extantStop){
		currentSimTime += simt->getTimeToNextVarEvent();
		keepSimulating = false;
	}
	
	return keepSimulating;
	
}


// =========================================================================
//	Simulator::processSimulation
// =========================================================================
void Simulator::processSimulation(bool scT){
	
	extantTree = new Tree(rando, extantStop);
	extantTree->setIncludeRootFossSampBool(includeRootFossils);
	simt->prepTreeForReconstruction();
	Node *simRoot = simt->getTreeRoot();
	extantTree->reconstructTreeFromSimulation(simRoot);
	if(scT)
		extantTree->scaleTreeDepthToValue(treescale);
	extantTree->setExtTreeExtRoot();
}

// =========================================================================
//	Simulator::processGSASimulation
// =========================================================================
void Simulator::processGSASimulation(){
	
	Tree *tt = new Tree(rando, gsaExtStop + simt->getNumExtinctLineages());
	tt->setIncludeRootFossSampBool(includeRootFossils);
	tt->setAllFossilSample(allFossilsSamp);
	simt->prepGSATreeForReconstruction();
	Node *simRoot = simt->getTreeRoot();
	tt->reconstructTreeFromSimulation(simRoot);
	gsaTreess.push_back(tt);

}


// =========================================================================
//	Simulator::getWholeTreeStringSim
// =========================================================================
string Simulator::getWholeTreeStringSim(){
	
	return simt->getWholeTreeString();
}

// =========================================================================
//	Simulator::getExtantTreeStringSim
// =========================================================================
string Simulator::getExtantTreeStringSim(){
	
	return extantTree->getWholeTreeString();
	
}

// =========================================================================
//	Simulator::getRawCalScaledTreeStringSim
// =========================================================================
string Simulator::getRawCalScaledTreeStringSim(){
	
	return extantTree->getCalScaleTreeDesc();
	
}


// =========================================================================
//	Simulator::getExtantTreeDepth
// =========================================================================
double Simulator::getExtantTreeDepth(){
	
	return extantTree->getTotalTreeDepth();
	
}



// =========================================================================
//	Simulator::setExtGSAStop
// =========================================================================
void Simulator::setExtGSAStop(int extStop, int gsaStop){
	
	extantStop = extStop;
	gsaExtStop = gsaStop;
}


// =========================================================================
//	Simulator::doPaleoSimStuff
// =========================================================================
void Simulator::doPaleoSimStuff(std::ofstream &dOut, double phiF, int numFossToSamp, bool rndFosSamp, double kpnu){
	simt->generatePaleoPriorFoss(dOut, phiF, numFossToSamp,rndFosSamp,kpnu);
}

// =========================================================================
//	Simulator::getPaleoSimNewickTreeDesc
// =========================================================================
string Simulator::getPaleoSimNewickTreeDesc(void){
	
	return simt->getPaleoPNewickDesc();
}

// =========================================================================
//	Simulator::getPaleoSimCalibrationDataString
// =========================================================================
string Simulator::getPaleoSimCalibrationDataString(void){
	
	return simt->getCalibrationDataString();
}



// =========================================================================
//	Simulator::getAllPaleoFossDistDataString
// =========================================================================
string Simulator::getAllPaleoFossDistDataString(void){
	
	return simt->getFossDistAncString();
}

// =========================================================================
//	Simulator::getTreeNumFossilExtantLinFossilString
// =========================================================================
string Simulator::getTreeNumFossilExtantLinFossilString(void){
	
	stringstream ss;
	int el = simt->getNumExLineageFossils();
	int tl = simt->getNumTotalFossilEvents();
	int nx = simt->getNumExtinctLineagesOnTree();
	int nl = simt->getNumExtPPTreeNodes() - 1;
	double ttl = simt->getTotalTreeLengthWholeT();
	double xtl = simt->getTotalLenghtOfExtinctLineages();
	ss << el << "\t" << tl << "\t" << (float)el / (float)tl << "\t";
	ss << nx << "\t" << nl << "\t" << (float)nx / (float)nl << "\t";
	ss << xtl << "\t" << ttl << "\t" << xtl / ttl << "\t";
	ss << simt->getExtantRootDepth() << "\n";
	return ss.str();
	
}



