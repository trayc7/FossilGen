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
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "engine.h"
#include "simulator.h"


using namespace std;



// =========================================================================
//	PhyEngine::PhyEngine
// =========================================================================
PhyEngine::PhyEngine(string of, int mt, double br, double dr, 
					 int nt, int ba, int da, double bp, double dp, 
					 double sc, double rp, int sd1, int sd2,
					 bool ppp, double frt, double sgal, double sgrt, bool doSG, 
					 int fsamp, bool inRt, bool rnsmp, double kpanu, bool allFos){
	outfilename = of;
	modelType = mt;
	extantStop = nt;
	stBirthRt = br;
	stDeathRt = dr;
	birthSampAlpha = ba;
	deathSampAlpha = da;
	bPriorAlpha = bp;
	dPriorAlpha = dp;
	numTrees = rp;
	treescale = sc;
	doPaleoPrior = ppp;
	fossRatePhi = frt;
	doScaleRanGamma = doSG;
	scaleGammaAlpha = sgal;
	scaleGammaRate = sgrt;
	numFossilSamples = fsamp;
	includeRootFossilSamp = inRt;
	rndFossSamp = rnsmp;
	kappaNu = kpanu;
	allFossils = allFos;
	if(sd1 > 0 && sd2 > 0)
		rando.setSeed(sd1, sd2);
	else
		rando.setSeed();
	seedType gs1, gs2;
	rando.getSeed(gs1, gs2);
	cout << "\nSeeds = {" << gs1 << ", " << gs2 << "}" << endl;
	if(treescale > 0.0)
		doScaleTree = true;
	else
		doScaleTree = false;
}

// =========================================================================
//	PhyEngine::~PhyEngine
// =========================================================================
PhyEngine::~PhyEngine(){
		
	
}





// =========================================================================
//	PhyEngine::doRunRun
// =========================================================================
void PhyEngine::doRunRun(){
	
	if(stBirthRt <= 0.0){
		cerr << "ERROR: The starting birth rate is set to " << stBirthRt << " --  This is not right!" << endl;
		cerr << "No trees were generated." << endl;
		exit(1);
	}
	doSimulatePaleoPriorFoss();
}	





// =========================================================================
//	PhyEngine::doSimulatePaleoPriorFoss
// =========================================================================

void PhyEngine::doSimulatePaleoPriorFoss(){
	
	ofstream datOut, fdOut, nfOut;
	string ofn = outfilename + ".fdist.dat";
	string ofnFD = outfilename + ".all_fdanc.dat";
	string numFFN = outfilename + ".nFoss.dat";
	datOut.open(ofn.c_str());
	fdOut.open(ofnFD.c_str());
	nfOut.open(numFFN.c_str());
	double sumTS = 0.0;
	for(int i=0; i<numTrees; i++){
		if(doScaleRanGamma){
			treescale = rando.gammaRv(scaleGammaAlpha, scaleGammaRate);
			sumTS += treescale;
		}
		Simulator *treesim = new Simulator(&rando, modelType, extantStop, treescale,
										   stBirthRt, stDeathRt, birthSampAlpha, deathSampAlpha,
										   bPriorAlpha, dPriorAlpha, includeRootFossilSamp,allFossils);
		treesim->setExtGSAStop(extantStop*100, extantStop);
		TreeInfo *ti = new TreeInfo(i);
		simtrees.push_back(ti);
		treesim->simulateTree();
		ti->setWholeTreeStringInfo(treesim->getWholeTreeStringSim());
		ti->setExtantTreeStringInfo(treesim->getExtantTreeStringSim());
		ti->setTimeTreeStringInfo(treesim->getRawCalScaledTreeStringSim());
		ti->setTreeDepthInfo(treesim->getExtantTreeDepth());
		treesim->doPaleoSimStuff(datOut,fossRatePhi,numFossilSamples,rndFossSamp,kappaNu);
		if(numFossilSamples > 0){
			ti->setPaleoFossFigTreeString(treesim->getPaleoSimNewickTreeDesc());
			ti->setPaleoFossCalibData(treesim->getPaleoSimCalibrationDataString());
			fdOut << treesim->getAllPaleoFossDistDataString();
		}
		nfOut << treesim->getTreeNumFossilExtantLinFossilString();
	}
	datOut.close();
	fdOut.close();
	nfOut.close();
	writeTreeFiles();
	calcAverageRootAges();
	if(doScaleRanGamma)
		cout << "  ** Average clock rate = " << (sumTS / (float)numTrees) << " **\n";
	
}



// =========================================================================
//	PhyEngine::writeTreeFiles
// =========================================================================
void PhyEngine::writeTreeFiles(){
	
	
	for(vector<TreeInfo *>::iterator p = simtrees.begin(); p != simtrees.end(); p++){
		(*p)->writeWholeTreeFileInfo(outfilename);
		(*p)->writeExtantTreeFileInfo(outfilename);
		(*p)->writeScaledDepthTreeFileInfo(outfilename);
		if(doPaleoPrior)
			(*p)->writePaleoFossilTreeFileInfo(outfilename);
	}
	
	
}

// =========================================================================
//	PhyEngine::findTreeByID
// =========================================================================
TreeInfo* PhyEngine::findTreeByID(int i){
	
	TreeInfo *q;
	for(vector<TreeInfo *>::iterator p = simtrees.begin(); p != simtrees.end(); p++){
		if((*p)->getTreeIDInfo() == i){
			q = (*p);
			break;
		}
	}
	return q;
}

// =========================================================================
//	PhyEngine::calcAverageRootAges
// =========================================================================
void PhyEngine::calcAverageRootAges(){
	
	ofstream out;
	out.open("Average_root_depths.out");
	double sumRH = 0.0;
	for(vector<TreeInfo *>::iterator p = simtrees.begin(); p != simtrees.end(); p++){
		sumRH += (*p)->getTreeDepthInfo();
		out << (*p)->getTreeDepthInfo() << "\n";
	}
	out.close();
	
	cout << "\n  ** Average root age = " << sumRH / numTrees << " **" << endl;
	
}



// ###################### TreeInfo ######################################################


// =========================================================================
//	TreeInfo::writeWholeTreeFileInfo
// =========================================================================
void TreeInfo::writeWholeTreeFileInfo(string ofp){
	
	string fn = ofp;
	ofstream out;
	stringstream tn;
	tn << tInd;
	
	fn += "_" + tn.str() + ".w.tre";
	out.open(fn.c_str());
	out << "#NEXUS\nbegin trees;\n    tree wholeT_" << tInd << " = ";
	out << wholeT << "\n";
	out << "end;";
	out.close();
}	


// =========================================================================
//	TreeInfo::writeExtantTreeFileInfo
// =========================================================================
void TreeInfo::writeExtantTreeFileInfo(string ofp){
	
	string fn = ofp;
	ofstream out;
	stringstream tn;
	tn << tInd;
	
	fn += "_" + tn.str() + ".ext.tre";
	out.open(fn.c_str());
	out << "#NEXUS\nbegin trees;\n    tree extantT_" << tInd << " = ";
	out << extantTree << "\n";
	out << "end;";
	out.close();
	fn = ofp + "_" + tn.str() + ".ext.phy";
	out.open(fn.c_str());
	out << extantTree << "\n";
	out.close();
	
}	

// =========================================================================
//	TreeInfo::writeScaledDepthTreeFileInfo
// =========================================================================
void TreeInfo::writeScaledDepthTreeFileInfo(string ofp){
	
	string fn = ofp;
	ofstream out;
	stringstream tn;
	tn << tInd;
	
	fn += "_" + tn.str() + ".SCL.tre";
	out.open(fn.c_str());
	out << "#NEXUS\nbegin trees;\n    tree timeTree" << tInd << " = ";
	out << timeScaleTree << "\n";
	out << "end;";
	out.close();
	//cout << timeScaleTree << endl;

}	

// =========================================================================
//	TreeInfo::writeFossilTreeFileInfo
// =========================================================================
void TreeInfo::writeFossilTreeFileInfo(string ofp){
	
	string fn = ofp;
	ofstream out;
	stringstream tn;
	tn << tInd;
	
	fn += "_" + tn.str() + ".WFOS.tre";
	out.open(fn.c_str());
	out << "#NEXUS\nbegin trees;\n    tree fossilTree" << tInd << " = ";
	out << fossilFig << "\n";
	out << "end;";
	out.close();
}	

// =========================================================================
//	TreeInfo::writeGSATreeFileInfo
// =========================================================================
void TreeInfo::writeGSATreeFileInfo(string ofp){
	
	string fn = ofp;
	ofstream out;
	stringstream tn;
	tn << tInd;
	
	fn += "_" + tn.str() + ".GSA.tre";
	out.open(fn.c_str());
	out << "#NEXUS\nbegin trees;\n";
	for(int t=0; t<gsaTree.size(); t++){
		out << "    tree timeTree" << tInd << "_" << t << " = ";
		out << gsaTree[t] << "\n";
	}

	out << "end;";
	out.close();
}	

// =========================================================================
//	TreeInfo::writePaleoFossilTreeFileInfo
// =========================================================================
void TreeInfo::writePaleoFossilTreeFileInfo(string ofp){
	
	string fn = ofp;
	ofstream out;
	stringstream tn;
	tn << tInd;
	
	fn += "_" + tn.str() + ".PPF.tre";
	out.open(fn.c_str());
	out << "#NEXUS\nbegin trees;\n    tree paleoFossilTree" << tInd << " = ";
	out << pfPaleoTreeFig << "\n";
	out << "end;";
	out.close();
	
	fn = ofp + "_" + tn.str() + ".ppf_info.cal";
	out.open(fn.c_str());
	out << pfCalibData << "\n";
	out.close();
}	


