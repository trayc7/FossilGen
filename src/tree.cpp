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
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include "tree.h"

using namespace std;

// =========================================================================
//	Fossil::Fossil
// =========================================================================
Fossil::Fossil(void){
	decN = NULL;
	ancN = NULL;
	age = -1.0;
	probSamp = -1.0;
	index = -1;
	isSampled = false;
	isCalib = false;
}


// end Fossil code
// =========================================================================




// =========================================================================
//	Node::Node
// =========================================================================
Node::Node(void){
	ldes = NULL;
	rdes = NULL;
	anc = NULL;
	sib = NULL;
	indx = -1;
	tname = "";
	flag = 0;
	isTip = false;
	branchLen = 0.0;
	nodeDepth = 0.0;
	birthRate = 0.0;
	deathRate = 0.0;
	birthTime = 0.0;
	deathTime = 0.0;
	isRoot = false;
	prSampFoss = 0.0;
	isIGFossil = false;
	isFoundFossil = false;
	isCalibratedNode = false;
	calibMinAge = 0.0;
	numDecTaxa = 0;
	randTestCals = NULL;
	
	fossPositions = NULL;
	fossSampProbs = NULL;
	numFossEvents = 0;
	numAllDecTaxa = 0;
	lDFD = -100.0;
	rDFD = -100.0;
	aFD = -100.0;
	nearFoss = -100.0;
	lDirF = NULL;
	rDirF = NULL;
	aDirF = NULL;
	lDFTime = -100.0;
	rDFTime = -100.0;
	aFTime = -100.0;
	nearTime = -1000.0;
	sampPaleoProb = 0.0;
	myFossils = NULL;
	if(!foundFossils.empty())
		foundFossils.clear();
	isExtRoot = false;
	calFossAge = -1.0;
	
	
}

// =========================================================================
//	Node::~Node
// =========================================================================
Node::~Node(){
	if(randTestCals != NULL)
		delete [] randTestCals;
	if(fossPositions != NULL)
		delete [] fossPositions;
	if(fossSampProbs != NULL)
		delete [] fossSampProbs;
	if(!foundFossils.empty())
		foundFossils.clear();
}

// =========================================================================
//	Node::nullNode
// =========================================================================
void Node::nullNode(void){
	ldes = NULL;
	rdes = NULL;
	anc = NULL;
	sib = NULL;
	indx = -1;
	tname = "";
	flag = 0;
	isTip = false;
	branchLen = 0.0;
	nodeDepth = 0.0;
	birthRate = 0.0;
	deathRate = 0.0;
	birthTime = 0.0;
	deathTime = 0.0;
	isRoot = false;
	prSampFoss = 0.0;
	isIGFossil = false;
	isFoundFossil = false;
	isCalibratedNode = false;
	calibMinAge = 0.0;
	numDecTaxa = 0;
	if(randTestCals != NULL)
		delete [] randTestCals;
	randTestCals = NULL;

	if(fossPositions != NULL)
		delete [] fossPositions;
	if(fossSampProbs != NULL)
		delete [] fossSampProbs;
	if(myFossils != NULL)
		delete [] myFossils;
	if(!foundFossils.empty())
		foundFossils.clear();
	
	fossPositions = NULL;
	fossSampProbs = NULL;
	numFossEvents = 0;
	numAllDecTaxa = 0;
	lDFD = -100.0;
	rDFD = -100.0;
	aFD = -100.0;
	nearFoss = -100.0;
	lDirF = NULL;
	rDirF = NULL;
	aDirF = NULL;
	lDFTime = -100.0;
	rDFTime = -100.0;
	aFTime = -100.0;
	nearTime = -1000.0;
	sampPaleoProb = 0.0;
	myFossils = NULL;
	isExtRoot = false;
	calFossAge = -1.0;

}

// =========================================================================
//	Node::nullPPFD
// =========================================================================
void Node::nullPPFD(void){
	if(fossPositions != NULL)
		delete [] fossPositions;
	if(fossSampProbs != NULL)
		delete [] fossSampProbs;
	if(myFossils != NULL)
		delete [] myFossils;
	if(!foundFossils.empty())
		foundFossils.clear();
	
	fossPositions = NULL;
	fossSampProbs = NULL;
	numFossEvents = 0;
	numAllDecTaxa = 0;
	lDFD = -100.0;
	rDFD = -100.0;
	aFD = -100.0;
	nearFoss = -100.0;
	lDirF = NULL;
	rDirF = NULL;
	aDirF = NULL;
	lDFTime = -100.0;
	rDFTime = -100.0;
	aFTime = -100.0;
	nearTime = -1000.0;
	sampPaleoProb = 0.0;
	isExtRoot = false;
	calFossAge = -1.0;

}


// =========================================================================
//	Node::initializeFossList
// =========================================================================
void Node::initializeFossList(int &ix, bool inclRoot){
	Node *extAnc = anc;
	int flag = extAnc->getFlag();
	while(flag < 2){
		extAnc = extAnc->getAnc();
		flag = extAnc->getFlag();
	}
	if(numFossEvents > 0){
		myFossils = new Fossil[numFossEvents];
		for(int i=0; i<numFossEvents; i++){
			Fossil *f = &myFossils[i];
			double fAge = getAgeOfFossilAtPosN(i);
			double fProb = fossSampProbs[i];
			f->setAge(fAge);
			f->setDecN(this);
			f->setTrueProbSamp(fProb);
			f->setFossIndex(ix);
			f->setAncN(extAnc);
			if(extAnc->getIsExtantRoot() && !inclRoot){
				f->setProbSamp(0.0);
			}
			else
				f->setProbSamp(fProb);
			ix++;
		}
		
	}
	else
		myFossils = NULL;
}

// =========================================================================
//	Node::getAgeForFossIndex
// =========================================================================
double Node::getAgeForFossIndex(int i){
	double rAge = -1.0;
	if(numFossEvents < 0){
		return rAge;
	}
	else{
		
		for(int i=0; i<numFossEvents; i++){
			Fossil *f = &myFossils[i];
			if(f->getFossIndex() == i)
				rAge = f->getAge();
		}
	}
	return rAge;
}
// =========================================================================
//	Node::getNumSampledMyFossils
// =========================================================================
int Node::getNumSampledMyFossils(void){
	
	int sum = 0;
	for(int i=0; i<numFossEvents; i++){
		Fossil *f = &myFossils[i];
		if(f->getIsFosSampled())
			sum++;
	}
	return sum;
}

// =========================================================================
//	Node::setSampFossAgesArray
// =========================================================================
void Node::setSampFossAgesArray(double *a){
	
	//int ns = getNumSampledMyFossils();
	int count = 0;
	for(int i=0; i<numFossEvents; i++){
		Fossil *f = &myFossils[i];
		if(f->getIsFosSampled()){
			a[count] = f->getAge();
			count++;
		}
	}	
}

// =========================================================================
//	Node::falseFossilCalibrating
// =========================================================================
void Node::falseFossilCalibrating(void){
	
	for(vector<Fossil *>::iterator f=foundFossils.begin(); f!=foundFossils.end(); f++){
		(*f)->setIsCalibrating(false);
	}
}



// end Node code 
// =========================================================================



// =========================================================================
//	Tree::Tree
// =========================================================================
Tree::Tree(MbRandom *p, double b0, double d0, int ext, double ct, int ba, int da,
		   double bp, double dp, bool inRt, bool allF){
	
	rando = p;
	allocated = 0;
	available = 0;
	capacity = 0;
	ntax = 0;
	numNodes = 0;
	numExtant = 0;
	numExtinct = 0;
	numTotalTips = 0;
	bRate0 = b0;
	dRate0 = d0;
	
	birthSampAlpha = ba;
	deathSampAlpha = da;
	bPriorAlpha = bp;
	dPriorAlpha = dp;
	nodes = NULL;
	rawTreeDepth = 1.0;
	scaleToVal = 1.0;
	ntax = ext;
	upPassSeq = NULL;
	initializeTree(ct);
	extUpPassSeq = NULL;
	extDownPassSeq = NULL;
	numExtPPTreeNodes = 0;
	numExtPPTreeTips = 0;
	numTotalFossilEvents = 0;
	paleoFossSampleNum = -1;
	allFDString = "";
	includeRootFoss = inRt;
	exLinEvents = 0;
	nExtinctLineages = 0;
	totalExtinctTL = 0.0;
	extantRootDepth = 0.0;
	sampAllFossils = allF;
	
}

// =========================================================================
//	Tree::Tree
// =========================================================================
Tree::Tree(MbRandom *p, int nt){
	rando = p;
	ntax = nt;
	numNodes = 2 * ntax - 1;
	nodes = new Node[numNodes];
	upPassSeq = new Node*[numNodes];
	for(int i=0; i<numNodes; i++){
		nodes[i].nullNode();
		nodes[i].setIndex(i);
	}
}


// =========================================================================
//	Tree::Tree
// =========================================================================
Tree::~Tree(){
	
	if(nodes != NULL){
		for(int i=0; i<numNodes; i++)
			nodes[i].nullNode();
		delete [] nodes;
	}
	if(upPassSeq != NULL)
		delete [] upPassSeq;
	if(extUpPassSeq != NULL)
		delete [] extUpPassSeq;
	if(extDownPassSeq != NULL)
		delete [] extDownPassSeq;
	if(!simulating.empty())
		simulating.clear();
}


// =========================================================================
//	Tree::initializeTree
// =========================================================================
void Tree::initializeTree(double curTime){
	
	currentTime = curTime;
	if(simulating.empty()){
		Node *p = new Node();
		root = p;
		simulating.push_back(p);
	}
	else{
		simulating.clear();
		Node *p = new Node();
		root = p;
		simulating.push_back(p);
	}
	
	root->setNBirthRate(bRate0);
	root->setNDeathRate(dRate0);
	root->setNodeBirthTime(curTime);
	root->setIsRoot(true);
	extants.push_back(root);
	numExtant = 1;
}

// =========================================================================
//	Tree::getTimeToNextERMEvent
// =========================================================================
double Tree::getTimeToNextERMEvent(){
	
	double sumrt = bRate0 + dRate0;
	double returnTime = 0.0;
	
	returnTime = -log(rando->uniformRv()) / (double(numExtant) * sumrt);
	return returnTime;
}


// =========================================================================
//	Tree::getTimeToNextVarEvent
// =========================================================================
double Tree::getTimeToNextVarEvent(){
	
	double sumrt = 0.0;
	double returnTime = 0.0;
	for(vector<Node *>::iterator p=extants.begin(); p!=extants.end(); p++){
		sumrt += (*p)->getNBirthRate();
		sumrt += (*p)->getNDeathRate();
	}
	returnTime = -log(rando->uniformRv()) / (double(numExtant) * sumrt);
	
	return returnTime;
}

// =========================================================================
//	Tree::ermEvent
// =========================================================================
void Tree::ermEvent(double ct){

	currentTime = ct;
	int nodeInd = rando->discreteUniformRv(0, numExtant-1);
	Node *eventN = (*(extants.begin() + nodeInd));
	double relBR = bRate0 / (bRate0 + dRate0);
	bool isBirth = ( rando->uniformRv() < relBR ? true : false );
	if(isBirth)
		lineageBirthEvenet(eventN, nodeInd);
	else
		lineageDeathEvent(eventN, nodeInd);
		
}	

// =========================================================================
//	Tree::varEvent
// =========================================================================
void Tree::varEvent(double ct){
	
	currentTime = ct;
	int numRates = 2 * numExtant;
	double *rates = new double[numRates];
	double *probs = new double[numRates];
	double sumrt = 0.0;
	int i = 0;
	for(vector<Node *>::iterator p=extants.begin(); p!=extants.end(); p++){
		rates[i] = (*p)->getNBirthRate();
		sumrt += rates[i];
		i++;
		rates[i] = (*p)->getNDeathRate();
		sumrt += rates[i];
		i++;
	}
	for(int j=0; j<numRates; j++){
		double pr = rates[j] / sumrt;
		if(j > 0)
			probs[j] = pr + probs[j-1];
		else
			probs[j] = pr;
	}
	
	if(probs[numRates-1] < 1.0){
		for(int j=0; j<numRates; j++)
			probs[j] /= probs[numRates-1];
	}
	double ran = rando->uniformRv();
	int outcome = 0;
	while(ran > probs[outcome])
		outcome++;
	int nodeInd = floor(outcome * 0.5);
	Node *eventN = (*(extants.begin()+nodeInd));
	if(outcome % 2)
		lineageDeathEvent(eventN, nodeInd);
	else
		lineageBirthEvenet(eventN, nodeInd);
	delete [] rates;
	delete [] probs;
}



// =========================================================================
//	Tree::setExtantPresentTime
// =========================================================================
void Tree::setExtantPresentTime(double ct){
	
	currentTime = ct;
	for(vector<Node *>::iterator p=extants.begin(); p!=extants.end(); p++){
		(*p)->setNodeDeathTime(ct);
		(*p)->setIsExtant(true);
	}
	setTreeBranchLengths();
	setTreeNodeDepths();
	setTreeTipNames();
	numNodes = (int)simulating.size();
}


// =========================================================================
//	Tree::lineageDeathEvent
// =========================================================================
void Tree::lineageDeathEvent(Node *p, int idx){
	
	p->setNodeDeathTime(currentTime);
	p->setIsExtant(false);
	p->setIsTip(true);
	extants.erase(extants.begin()+idx);
	numExtant = (int)extants.size();
	numExtinct += 1;
}

// =========================================================================
//	Tree::lineageBirthEvenet
// =========================================================================
void Tree::lineageBirthEvenet(Node *p, int idx){
	
	Node *lefty, *sis;
	lefty = new Node();
	sis = new Node();
	extants.erase(extants.begin()+idx);
	setNewLineageInfo(p, lefty, sis);
	mutateLineage(p, lefty);
	mutateLineage(p, sis);
}


// =========================================================================
//	Tree::setNewLineageInfo
// =========================================================================
void Tree::setNewLineageInfo(Node *par, Node *offsp, Node *sib){
	
	par->setLdes(offsp);
	par->setRdes(sib);
	par->setNodeDeathTime(currentTime);
	par->setIsTip(false);
	par->setIsExtant(false);
	
	offsp->setLdes(NULL);
	offsp->setRdes(NULL);
	offsp->setSib(sib);
	offsp->setAnc(par);
	offsp->setNodeBirthTime(currentTime);
	offsp->setIsTip(true);
	offsp->setIsExtant(true);
	
	sib->setLdes(NULL);
	sib->setRdes(NULL);
	sib->setSib(offsp);
	sib->setAnc(par);
	sib->setNodeBirthTime(currentTime);
	sib->setIsTip(true);
	sib->setIsExtant(true);
	
	extants.push_back(offsp);
	extants.push_back(sib);
	simulating.push_back(offsp);
	simulating.push_back(sib);
	numExtant = (int)extants.size();
}


// =========================================================================
//	Tree::mutateLineage
// =========================================================================
void Tree::mutateLineage(Node *par, Node *dec){
	
	double pBRate = par->getNBirthRate();
	double pDRate = par->getNDeathRate();
	if(modelsim == 2)
		gammaMutateLineage(dec, pBRate, pDRate);
	else{
		dec->setNBirthRate(pBRate);
		dec->setNDeathRate(pDRate);
	}
}


// =========================================================================
//	Tree::gammaMutateLineage
// =========================================================================
void Tree::gammaMutateLineage(Node *offsp, double pb, double pd){
	
	double r = 0.0;
	double accept = 0.0;
	double newRate;
	double grand;
	while(r >= accept){
		grand = rando->gammaRv(birthSampAlpha, birthSampAlpha);
		newRate = grand * pb;
		r = rando->uniformRv();
		accept = exp((bPriorAlpha-1.0)*(log(newRate)-log(pb))-(bPriorAlpha*newRate)+(bPriorAlpha*pb));
	}
	offsp->setNBirthRate(newRate);
	r = 0.0;
	accept = 0.0;
	while(r >= accept){
		grand = rando->gammaRv(deathSampAlpha, deathSampAlpha);
		newRate = grand * pd;
		r = rando->uniformRv();
		accept = exp((dPriorAlpha-1.0)*(log(newRate)-log(pd))-(dPriorAlpha*newRate)+(dPriorAlpha*pd));
	}
	offsp->setNDeathRate(newRate);
	
}


// =========================================================================
//	Tree::zeroAllFlags
// =========================================================================
void Tree::zeroAllFlags(){
	
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++)
		(*p)->setFlag(0);
}


// =========================================================================
//	Tree::setWholeTreeFlags
// =========================================================================
void Tree::setWholeTreeFlags(){
	
	zeroAllFlags();
	numTotalTips = 0;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		if((*p)->getIsTip()){
			(*p)->setFlag(1);
			numTotalTips++;
		}
	}
	setSampleFromFlags();
}

// =========================================================================
//	Tree::setExtantTreeFlags
// =========================================================================
void Tree::setExtantTreeFlags(){
	
	zeroAllFlags();
	for(vector<Node *>::iterator p=extants.begin(); p!=extants.end(); p++)
		(*p)->setFlag(1);
	setSampleFromFlags();
}

// =========================================================================
//	Tree::setGSATipTreeFlags
// =========================================================================
void Tree::setGSATipTreeFlags(){
	
	zeroAllFlags();
	numTotalTips = 0;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		if((*p)->getIsTip()){
			(*p)->setFlag(1);
			numTotalTips++;
		}
	}
	setSampleFromFlags();
	
}

// =========================================================================
//	Tree::popExtantVec
// =========================================================================
void Tree::popExtantVec(){
	
	//cout << simulating.size() << endl;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		simulating.push_back(p);
		if(p->getIsTip()){
			if(p->getIsExtant()){
				extants.push_back(p);
			}
		}
	}
	
}




// =========================================================================
//	Tree::setSampleFromFlags
// =========================================================================
void Tree::setSampleFromFlags(){

	
	int flag;
	Node *q = NULL;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		if((*p)->getIsTip()){
			flag = (*p)->getFlag();
			q = (*p);
			if(flag == 1){
				do{
					q = q->getAnc();
					flag = q->getFlag();
					flag++;
					q->setFlag(flag);
				}while (q->getIsRoot() == false && flag < 2);
			}
		}
	}
}

	
// =========================================================================
//	Tree::setWholeTreeExtantRoot 
// =========================================================================
// putting this after creating the extant tree
void Tree::setWholeTreeExtantRoot(){
	setExtantTreeFlags();
	for(int i=0; i<numNodes; i++){
		Node *p = upPassSeq[i];
		if(p->getFlag() == 2){
			extantRoot = p;
			extantRoot->setIsExtantRoot(true);
			break;
		}
		
	}
}




// =========================================================================
//	Tree::setTreeBranchLengths
// =========================================================================
void Tree::setTreeBranchLengths(){
	
	double btime;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		btime = (*p)->getNodeDeathTime() - (*p)->getNodeBirthTime();
		(*p)->setBranchLen(btime);
	}
}

// =========================================================================
//	Tree::setTreeNodeDepths
// =========================================================================
void Tree::setTreeNodeDepths(){
	
	double myd;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		myd = currentTime - (*p)->getNodeDeathTime();
		(*p)->setNodeDepth(myd);
	}
}


// =========================================================================
//	Tree::setTreeTipNames
// =========================================================================
void Tree::setTreeTipNames(){
	
	int extantIt = 0;
	int tipI = ntax;
	bool recNames = true;
	if(recNames){
		setRecTreeTipNames(root, extantIt, tipI);
	}
	else {
		for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
			stringstream tn;
			if((*p)->getIsTip()){
				if((*p)->getIsExtant()){
					tn << extantIt + 1;
					string name = "T" + tn.str();
					(*p)->setTipName(name);
					extantIt++;
				}
				else{
					tn << tipI + 1;
					string name = "X" + tn.str();
					(*p)->setTipName(name);
					tipI++;
				}
			}
		}
	}
}


// =========================================================================
//	Tree::setRecTreeTipNames
// =========================================================================
void Tree::setRecTreeTipNames(Node *inN, int &extIt, int &tipIt){
	
	if (inN != NULL){
		if(inN->getIsTip()){
			stringstream tn;
			if(inN->getIsExtant()){
				tn << extIt + 1;
				string name = "T" + tn.str();
				inN->setTipName(name);
				extIt++;
			}
			else{
				tn << tipIt + 1;
				string name = "X" + tn.str();
				inN->setTipName(name);
				tipIt++;
			}
		}
		else{
			setRecTreeTipNames(inN->getLdes(), extIt, tipIt);
			setRecTreeTipNames(inN->getRdes(), extIt, tipIt);
		}
	}
}


// =========================================================================
//	Tree::getWholeTreeString
// =========================================================================
string Tree::getWholeTreeString(){
	
	return getIncWholeTreeNewick();
}

// =========================================================================
//	Tree::prepTreeForReconstruction
// =========================================================================
void Tree::prepTreeForReconstruction(){
	
	setExtantTreeFlags();
}

// =========================================================================
//	Tree::prepGSATreeForReconstruction
// =========================================================================
void Tree::prepGSATreeForReconstruction(){
	
	setGSATipTreeFlags();
}


// =========================================================================
//	Tree::scaleTreeDepthToValue
// =========================================================================
void Tree::scaleTreeDepthToValue(double scVal){
	
	double depth = 0.0;
	int k = 0;
	Node *p = &nodes[k];
	while(!p->getIsTip())
		p = &nodes[k++];
	while(p != root){
		depth += p->getBranchLen();
		p = p->getAnc();
	}
	rawTreeDepth = depth;
	scaleToVal = scVal;
	double scaler = scVal / depth;

	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsRoot() == false){
			double bl = p->getBranchLen();
			p->setBranchLen(bl * scaler);
		}
	}
}


// =========================================================================
//	Tree::reconstructTreeFromSimulation
// =========================================================================
void Tree::reconstructTreeFromSimulation(Node *oRoot){
	
	Node *p = NULL;
	int nextIntNode = ntax;
	int tipc = 0;
	reconstructLineageFromSimulation(p, oRoot, nextIntNode, tipc);
	getRevDownPassSequence();
}

// =========================================================================
//	Tree::reconstructLineageFromSimulation
// =========================================================================
int Tree::reconstructLineageFromSimulation(Node *myN, Node *oN, int &nextIntN, int &tipcounter){
	
	Node *p = NULL;
	double brlen = oN->getBranchLen();
	double depth = oN->getNodeDepth();
	int oFlag = oN->getFlag();
	bool rootN = oN->getIsRoot();
	if(oN->getIsTip() && oFlag == 1){
		Node *oAnc = oN->getAnc();
		int aF = oAnc->getFlag();
		if(aF == 1){
			brlen += oAnc->getBranchLen();
			while(!oAnc->getIsRoot() && aF < 2){
				oAnc = oAnc->getAnc();
				aF = oAnc->getFlag();
				if(aF == 1)
					brlen += oAnc->getBranchLen();
			}
		}
		p = &nodes[tipcounter++];
		p->setIsTip(true);
		p->setIsExtant(oN->getIsExtant());
		if(!p->getIsExtant()) p->setIsExtant(false);
		p->setTipName(oN->getTipName());
		p->setBranchLen(brlen);
		p->setNodeDepth(depth);
		p->setNBirthRate(oN->getNBirthRate());
		p->setNDeathRate(oN->getNDeathRate());
		p->setPrSampFossil(oN->getPrSampFossil());
		p->setIsFoundFossil(oN->getIsFoundFossil());
		p->setCalibrationMinAge(oN->getCalibrationMinAge());
		p->setAnc(myN);
		if(myN->getLdes() == NULL)
			myN->setLdes(p);
		else if(myN->getRdes() == NULL)
			myN->setRdes(p);
		else{
			cerr << "ERROR: Problem adding a tip to the tree" << endl;
			exit(1);
		}
		return tipcounter;
	}
	
	else{
		if(oFlag > 1){
			Node *s1 = &nodes[nextIntN++];
			if(oN->getLdes()->getFlag() > 0)
				reconstructLineageFromSimulation(s1, oN->getLdes(), nextIntN, tipcounter);
			if(oN->getRdes()->getFlag() > 0)
				reconstructLineageFromSimulation(s1, oN->getRdes(), nextIntN, tipcounter);
			if(rootN == false){
				Node *oAnc = oN->getAnc();
				int aF = oAnc->getFlag();
				if(aF == 1){
					brlen += oAnc->getBranchLen();
					while(!oAnc->getIsRoot() && aF < 2){
						oAnc = oAnc->getAnc();
						aF = oAnc->getFlag();
						if(aF == 1)
							brlen += oAnc->getBranchLen();
					}
				}
				if(myN != NULL){
					s1->setBranchLen(brlen);
					s1->setNodeDepth(depth);
					s1->setNBirthRate(oN->getNBirthRate());
					s1->setNDeathRate(oN->getNDeathRate());
					s1->setPrSampFossil(oN->getPrSampFossil());
					s1->setCalibrationMinAge(oN->getCalibrationMinAge());
					s1->setIsCalibratedNode(oN->getIsCalibratedNode());
					s1->setAnc(myN);
					if(myN->getLdes() == NULL)
						myN->setLdes(s1);
					else if(myN->getRdes() == NULL)
						myN->setRdes(s1);
					else{
						cerr << "ERROR: Problem adding a internal to the tree" << endl;
						exit(1);
					}
				}
				else{
					s1->setIsRoot(true);
					root = s1;
					s1->setBranchLen(brlen);
					s1->setNodeDepth(depth);
					s1->setNBirthRate(oN->getNBirthRate());
					s1->setNDeathRate(oN->getNDeathRate());
					s1->setPrSampFossil(oN->getPrSampFossil());
					s1->setCalibrationMinAge(oN->getCalibrationMinAge());
					s1->setIsCalibratedNode(oN->getIsCalibratedNode());
				}
			}
			else if(rootN == true){
				s1->setIsRoot(true);
				root = s1;
				s1->setBranchLen(0.0);
				s1->setNodeDepth(depth);
				s1->setNBirthRate(oN->getNBirthRate());
				s1->setNDeathRate(oN->getNDeathRate());
				s1->setCalibrationMinAge(oN->getCalibrationMinAge());
				s1->setIsCalibratedNode(oN->getIsCalibratedNode());
			}
		}
		
		else if(oFlag == 1){
			if(oN->getRdes()->getFlag() == 0 && oN->getLdes()->getFlag() > 0)
				reconstructLineageFromSimulation(myN, oN->getLdes(), nextIntN, tipcounter);
			else
				reconstructLineageFromSimulation(myN, oN->getRdes(), nextIntN, tipcounter);
		}
	}
	return tipcounter;
}	


// =========================================================================
//	Tree::passUp
// =========================================================================
void Tree::passUp(Node *p, int *x) {
	
	if(p != NULL){
		upPassSeq[(*x)++] = p;
		passUp(p->getLdes(), x);
		passUp(p->getRdes(), x);
	}
}



// =========================================================================
//	Tree::getRevDownPassSequence
// =========================================================================
void Tree::getRevDownPassSequence(void) {
	
	int x = 0;
	passUp(root, &x);
}


// =========================================================================
//	Tree::getNewickTreeDesc
// =========================================================================
string Tree::getNewickTreeDesc(void) {
	
	stringstream ss;
	writeTreeLineage(root, ss);
	ss << ";";
	string treeDesc = ss.str();
	return treeDesc;
}

// =========================================================================
//	Tree::getIncWholeTreeNewick
// =========================================================================
string Tree::getIncWholeTreeNewick(void) {
	
	stringstream ss;
	writeTreeLineage(extantRoot, ss);
	ss << ";";
	string treeDesc = ss.str();
	return treeDesc;
}

// =========================================================================
//	Tree::writeTreeLineage
// =========================================================================
void Tree::writeTreeLineage(Node *p, stringstream &ss) {
	
	if (p != NULL){
		if(p->getLdes() == NULL) // tip
			ss << p->getTipName();
		else{
			ss << "(";
			writeTreeLineage(p->getLdes(), ss);
			ss << ":" << p->getLdes()->getBranchLen();
			ss << ",";
			writeTreeLineage(p->getRdes(), ss);
			ss << ":" << p->getRdes()->getBranchLen();
			ss << ")";
		}
	}
}



// =========================================================================
//	Tree::getCalScaleTreeDesc
// =========================================================================
string Tree::getCalScaleTreeDesc(void) {
	
	stringstream ss;
	writeCalScaleTreeLineage(root, ss);
	ss << ";";
	string treeDesc = ss.str();
	return treeDesc;
}

// =========================================================================
//	Tree::writeCalScaleTreeLineage
// =========================================================================
void Tree::writeCalScaleTreeLineage(Node *p, stringstream &ss) {
	
	if (p != NULL){
		if(p->getLdes() == NULL) // tip
			ss << p->getTipName();
		else{
			ss << "(";
			writeCalScaleTreeLineage(p->getLdes(), ss);
			ss << ":" << (p->getLdes()->getBranchLen() / scaleToVal)*rawTreeDepth;
			ss << ",";
			writeCalScaleTreeLineage(p->getRdes(), ss);
			ss << ":" << (p->getRdes()->getBranchLen() / scaleToVal)*rawTreeDepth;
			ss << ")";
		}
	}
}

// =========================================================================
//	Tree::findLeftLeafOfNode
// =========================================================================
string Tree::findLeftLeafOfNode(Node *p) {
	
	Node *d = p->getLdes();
	while(d->getIsTip() == false)
		d = d->getLdes();
	return d->getTipName();
}	

// =========================================================================
//	Tree::findRightLeafOfNode
// =========================================================================
string Tree::findRightLeafOfNode(Node *p) {
	
	Node *d = p->getRdes();
	while(d->getIsTip() == false)
		d = d->getLdes();
	return d->getTipName();
}	

// =========================================================================
//	Tree::getExtantTreeRoot
// =========================================================================
Node* Tree::getExtantTreeRoot(Node *p) {
	
	if(p->getIsTip())
		return NULL;
	else{
		int flag = p->getFlag();
		if(flag == 2)
			return p;
		else{
			Node *q = getExtantTreeRoot(p->getLdes());
			if(q == NULL)
				q = getExtantTreeRoot(p->getRdes());
			return q;
		}
	}
	return NULL;
}	




// =========================================================================
//	Tree::getFossilNewickTreeDesc
// =========================================================================
string Tree::getFossilNewickTreeDesc(void) {
	
	stringstream ss;
	writeFossilTreeLineage(root, ss);
	ss << ";";
	string treeDesc = ss.str();
	return treeDesc;
}

// =========================================================================
//	Tree::writeFossilTreeLineage
// =========================================================================
void Tree::writeFossilTreeLineage(Node *p, stringstream &ss) {
	
	if (p != NULL){
		if(p->getLdes() == NULL) // tip
			ss << p->getTipName();
		else{
			ss << "(";
			writeFossilTreeLineage(p->getLdes(), ss);
			ss << "[&fossil=" << p->getLdes()->getPrSampFossil() << "]";
			ss << ":" << p->getLdes()->getBranchLen();
			ss << ",";
			writeFossilTreeLineage(p->getRdes(), ss);
			ss << "[&fossil=" << p->getRdes()->getPrSampFossil() << "]";
			ss << ":" << p->getRdes()->getBranchLen();
			ss << ")";
		}
	}
}



// =========================================================================
//	Tree::setFossilExtantTreeFlags
// =========================================================================
int Tree::setFossilExtantTreeFlags(){
	
	zeroAllFlags();
	int nTaxInT = 0;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		if((*p)->getIsTip()){
			if((*p)->getIsExtant()){
				(*p)->setFlag(1);
				nTaxInT++;
			}
			else if((*p)->getIsFoundFossil()){
				(*p)->setFlag(1);
				nTaxInT++;
			}
			else
				(*p)->setFlag(0);
		}
	}
	setSampleFromFlags();
	return nTaxInT;
}

// =========================================================================
//	Tree::getTreeNumDecTaxa
// =========================================================================
int Tree::getTreeNumDecTaxa(){
	return getNodeNumDecTaxa(root);
}
	

// =========================================================================
//	Tree::getNodeNumDecTaxa
// =========================================================================
int Tree::getNodeNumDecTaxa(Node *p){
	
	int mydecs = 0;
	if(p->getIsTip()){
		if(p->getIsExtant())
			return 1;
	}
	else{
		mydecs = getNodeNumDecTaxa(p->getLdes());
		mydecs += getNodeNumDecTaxa(p->getRdes());
		p->setNumDecExtantTax(mydecs);
		return mydecs;
	}
	return mydecs;
}	


// =========================================================================
//	Tree::setAllTreeNodeIDs
// =========================================================================
void Tree::setAllTreeNodeIDs(){
	
	int tipsIDs = 0;
	int nodeIDs = numTotalTips;
	setTreeNodeIDs(root, tipsIDs, nodeIDs);
}	

// =========================================================================
//	Tree::setTreeNodeIDs
// =========================================================================
void Tree::setTreeNodeIDs(Node *p, int &tipID, int &nID){
	
	if(p->getIsTip()){
		p->setIndex(tipID);
		tipID++;
	}
	else{
		p->setIndex(nID);
		nID++;
		setTreeNodeIDs(p->getLdes(), tipID, nID);
		setTreeNodeIDs(p->getRdes(), tipID, nID);
	}
}



// =========================================================================
//	Tree::setExtantFlagsInFossilTree
// =========================================================================
void Tree::setExtantFlagsInFossilTree(){
	
	for(int i=0; i<numNodes; i++)
		nodes[i].setFlag(0);
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsTip() && p->getIsExtant())
			p->setFlag(1);
	}
	int flag = 0;
	Node *q;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsTip()){
			flag = p->getFlag();
			q = p;
			if(flag == 1){
				do{
					q = q->getAnc();
					flag = q->getFlag();
					flag++;
					q->setFlag(flag);
				}while(q->getIsRoot() == false && flag < 2);
			}
		}
	}
}

// =========================================================================
//	Tree::findLeftLeafOfFossilNode
// =========================================================================
string Tree::findLeftLeafOfFossilNode(Node *p) {
	
	Node *d = p->getLdes();
	while(d->getIsTip() == false){
		int flag = d->getLdes()->getFlag();
		if(flag == 0)
			d = d->getRdes();
		else
			d = d->getLdes();
	}
	return d->getTipName();
}	

// =========================================================================
//	Tree::findRightLeafOfFossilNode
// =========================================================================
string Tree::findRightLeafOfFossilNode(Node *p) {
	
	Node *d = p->getRdes();
	while(d->getIsTip() == false){
		int flag = d->getLdes()->getFlag();
		if(flag == 0)
			d = d->getRdes();
		else
			d = d->getLdes();
	}
	return d->getTipName();
}	

// =========================================================================
//	Tree::getNumberOfCalibNodes
// =========================================================================
int Tree::getNumberOfCalibNodes(){

	int ncal = 0;
	for(vector<Node *>::iterator p=simulating.begin(); p!=simulating.end(); p++){
		if((*p)->getIsCalibratedNode())
			ncal++;
	}
	return ncal;
}




// =========================================================================
//	Tree::getNodeTimeLowerBound
// =========================================================================
double Tree::getNodeTimeLowerBound(Node *p){
	
	double lb = p->getLdes()->getNodeDepth();
	if(p->getRdes()->getNodeDepth() > lb)
		lb = p->getRdes()->getNodeDepth();
	return lb;
}


// =========================================================================
//	Tree::quickSort
// =========================================================================

void Tree::quickSort(double numbers[], int arrSize){
	
	qSort(numbers, 0, arrSize-1);
	
}




// =========================================================================
//	Tree::qSort
// =========================================================================

void Tree::qSort(double numbers[], int left, int right){
	
	double pivot;
	int lHold, rHold;
	
	lHold = left;
	rHold = right;
	
	pivot = numbers[left];
	
	while(left < right){
		
		while((numbers[right] >= pivot) && (left < right))
			right--;
		if(left != right){
			
			numbers[left] = numbers[right];
			left++;
		}
		
		while((numbers[left] <= pivot) && (left < right))
			left++;
		
		if(left != right){
			numbers[right] = numbers[left];
			right--;
		}
	}
	
	numbers[left] = pivot;
	pivot = left;
	left = lHold;
	right = rHold;
	
	if(left < pivot) qSort(numbers, left, pivot - 1);
	
	if(right > pivot) qSort(numbers, pivot + 1, right);
}



// =========================================================================
//	Tree::generatePaleoPriorFoss
// =========================================================================

void Tree::generatePaleoPriorFoss(ofstream &dOut, double phi, int sampNum, bool rfs, double kpnu){
	
	rndFossSamp = rfs;
	kappaNuSamp = kpnu;
	paleoFossSampleNum = sampNum;
	if(paleoFossSampleNum < 2) 
		paleoFossSampleNum = 2;
	double lambda = phi;
	nullNodePPFossInfo();
	numExtPPTreeTips = countNTipsInExtFosTree();
	numExtPPTreeNodes = numExtPPTreeTips * 2 - 1;
	getRevDPSeqExtRootTree();
	getDownPSeqExtRootTree();
	setExtantFlagsInFossilTree();
	totETL = getTotalTreeLength();
	
	int nf = 0;
	extantRoot->setIsExtantRoot(true);
	for(int i=0; i<numExtPPTreeNodes; i++){ 
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			double nBrLen = p->getBranchLen();
			int numEvents = rando->poissonRv(nBrLen * lambda);
			
			if(numEvents > 0){
				nf += numEvents;
				double *tempFossPos = new double[numEvents];
				for(int j=0; j < numEvents; j++){
					double randPos = rando->uniformRv(0, nBrLen);
					tempFossPos[j] = randPos;
				}
				quickSort(tempFossPos, numEvents);
				
				p->setNumPPFossEvents(numEvents);
				p->initializeFossPos(numEvents);
				for(int j=0; j < numEvents; j++){
					p->setFossPosN(j, tempFossPos[j]);
				}
				delete [] tempFossPos;
			}
		}
	}
	
	
	

	getNodeNearestPPFoss(extantRoot, NULL);
	setRecNodeFosPositions();
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		int flag = p->getFlag();
		if(flag == 2){
			double myD = p->getNearFDist();
			if(myD > 0)
				dOut << myD << "\n";
		}
	}
	
	
	if(sampAllFossils){
		paleoFossSampleNum = nf;
		evolveFossSamplingProbArr();
	}
	else if(sampNum > 0){
		if(rndFossSamp) 
			evolveFossRandomSamplingProbArr();
		else
			evolveFossSamplingProbArr();
	}
	calcNumAncDecFossilsExtantTree();

}




// =========================================================================
//	Tree::getNodeNearestPPFoss
// =========================================================================

Node* Tree::getNodeNearestPPFoss(Node *p, Node *a){
	
	if(p != NULL){
		if(p->getIsTip()){
			int nEvOnBr = p->getNumPPFossEvents();
			p->setLDirFoss(NULL);
			p->setRDirFoss(NULL);
			if(nEvOnBr > 0){
				p->setADirFoss(p);
				return p;
			}
			else{
				p->setADirFoss(a);
				return NULL;
			}
		}
		else{
			int nEvOnBr = p->getNumPPFossEvents();
			if(nEvOnBr > 0){
				Node *myLDF = getNodeNearestPPFoss(p->getLdes(), p);
				p->setLDirFoss(myLDF);
				Node *myRDF = getNodeNearestPPFoss(p->getRdes(), p);
				p->setRDirFoss(myRDF);
				p->setADirFoss(p);
				return p;
			}
			else{
				
				Node *myLDF = getNodeNearestPPFoss(p->getLdes(), a);
				p->setLDirFoss(myLDF);
				Node *myRDF = getNodeNearestPPFoss(p->getRdes(), a);
				p->setRDirFoss(myRDF);
				p->setADirFoss(a);
				double lFAge = 0.0;
				double rFAge = 0.0;
				if(myLDF != NULL)
					lFAge = myLDF->getAgeOfFossilAtPosN(0);
				if(myRDF != NULL)
					rFAge = myRDF->getAgeOfFossilAtPosN(0);
				

				if(lFAge > 0.0 || rFAge > 0.0){
					if(lFAge > rFAge) 
						return myLDF;
					else if(rFAge > lFAge)
						return myRDF; 
					else
						return NULL;
				}
				else return NULL;
			}
		}
		
	}
	return NULL;

	
}

// =========================================================================
//	Tree::setRecNodeFosPositions
// =========================================================================
void Tree::setRecNodeFosPositions(void) {
	
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		Node *lF = p->getLDirFoss();
		Node *rF = p->getRDirFoss();
		Node *aF = p->getADirFoss();
		double myND = p->getNodeDepth();
		if(lF != NULL){
			double lFosA = lF->getAgeOfFossilAtPosN(0);
			double dist = myND - lFosA; 
			p->setLDFDist(dist);
			p->setLDFAge(lFosA);
			
		}
		else{ 
			p->setLDFDist(-200.00);
			p->setLDFAge(-200.00);
		}
		if(rF != NULL){
			double rFosA = rF->getAgeOfFossilAtPosN(0);
			double dist = myND - rFosA; 
			p->setRDFDist(dist);
			p->setRDFAge(rFosA);
		}
		else{ 
			p->setRDFDist(-200.00);
			p->setRDFAge(-200.00);
		}
		if(aF != NULL){
			int ancNFoss = aF->getNumPPFossEvents();
			if(ancNFoss > 0){
				double ancFosTime = aF->getAgeOfFossilAtPosN(ancNFoss - 1);
				double dist = ancFosTime - myND;
				p->setAFDist(dist);
				p->setAFAge(ancFosTime);
			}
			else
				cout << "something bad" << endl;
		}
		else{ 
			p->setAFDist(-200.00);
			p->setAFAge(-200.00);
		}
		getNodesDistToClosestFoss(p);
	}
}

// =========================================================================
//	Tree::getDistBWTwoNodes
// =========================================================================
double Tree::getDistBWTwoNodes(Node *p, Node *q, int updown) {
	double sum = 0.0;
	if(updown == 0){ 
		Node *tmp = p->getAnc();
		sum = p->getBranchLen();
		while(tmp != q){
			sum += tmp->getBranchLen();
			tmp = tmp->getAnc();
		}
		int qNFoss = q->getNumPPFossEvents();
		if(qNFoss > 0){
			double qFD = q->getBranchLen() - q->getFossPosN(qNFoss - 1);
			sum += qFD;
		}
		else{
			cout << "something wrong" << endl;
		}
		
	}
	else if(updown > 0){ 
		Node *tmp = q->getAnc();
		sum = 0.0;
		int qNFoss = q->getNumPPFossEvents();
		
		if(qNFoss > 0){
			sum = q->getFossPosN(0);
		}
		else{
			cout << "something wrong" << endl;
		}
		while(tmp != p){
			sum += tmp->getBranchLen();
			tmp = tmp->getAnc();
		}
		
	}
	return sum;
}




// =========================================================================
//	Tree::getNodesDistToClosestFoss
// =========================================================================
double Tree::getNodesDistToClosestFoss(Node *p) {
	double aD = p->getAFDist();
	double lD = p->getLDFDist();
	double rD = p->getRDFDist();
	double ttl = totETL;
	
	double aG = p->getAFAge();
	double lG = p->getLDFAge();
	double rG = p->getRDFAge();
	
	double closestA = -1000.0;
	//cout << ttl << endl;
	double lowest = ttl;
	if(aD > 0.0 && aD < lowest){
		lowest = aD;
		closestA = aG;
	}
	if(lD > 0.0 && lD < lowest){
		lowest = lD;
		closestA = lG;
	}
	if(rD > 0.0 && rD < lowest){
		lowest = rD;
		closestA = rG;
	}
	if(lowest == ttl){
		lowest = -200.00;
		closestA = -5000.0;
	}
	p->setNearFDist(lowest);
	p->setNearestFAge(closestA);
	return lowest;
}



// =========================================================================
//	Tree::getRevDPSeqExtRootTree
// =========================================================================
void Tree::getRevDPSeqExtRootTree(void) {
	if(extUpPassSeq != NULL)
		delete [] extUpPassSeq;
	extUpPassSeq = new Node*[numExtPPTreeNodes];
	int x = 0;
	passExtTUp(extantRoot, &x);
}

// =========================================================================
//	Tree::passExtTUp
// =========================================================================
void Tree::passExtTUp(Node *p, int *x) {
	
	
	if(p != NULL){
		extUpPassSeq[(*x)++] = p;
		passExtTUp(p->getLdes(), x);
		passExtTUp(p->getRdes(), x);
	}
}

// =========================================================================
//	Tree::getDownPSeqExtRootTree
// =========================================================================
void Tree::getDownPSeqExtRootTree(void) {
	if(extDownPassSeq != NULL)
		delete [] extDownPassSeq;
	extDownPassSeq = new Node*[numExtPPTreeNodes];
	int x = 0;
	passExtTDown(extantRoot, &x);
}

// =========================================================================
//	Tree::passExtTDown
// =========================================================================
void Tree::passExtTDown(Node *p, int *x) {
	
	
	if(p != NULL){
		passExtTDown(p->getLdes(), x);
		passExtTDown(p->getRdes(), x);
		extDownPassSeq[(*x)++] = p;
	}
}

// =========================================================================
//	Tree::countNTipsInExtFosTree
// =========================================================================
int Tree::countNTipsInExtFosTree(void){
	return setNumTipDec(extantRoot);
}


// =========================================================================
//	Tree::setNumTipDec
// =========================================================================
int Tree::setNumTipDec(Node *p){
	
	int mydecs = 0;
	if(p->getIsTip()){
		return 1;
	}
	else{
		mydecs = setNumTipDec(p->getLdes());
		mydecs += setNumTipDec(p->getRdes());
		p->setNumAllDecExtantTax(mydecs);
		return mydecs;
	}
	return mydecs;
}	



// =========================================================================
//	Tree::getTotalTreeLength
// =========================================================================
double Tree::getTotalTreeLength(void){
	double sum = 0.0;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		sum += p->getBranchLen();
	}
	return sum;
	
}

// =========================================================================
//	Tree::nullNodePPFossInfo
// =========================================================================
void Tree::nullNodePPFossInfo(void){
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		p->nullPPFD();
	}
	
}

// =========================================================================
//	Tree::getPaleoPNewickDesc
// =========================================================================
string Tree::getPaleoPNewickDesc(void) {
	
	stringstream ss;
	writePaleoPLineage(extantRoot, ss);
	ss << ";";
	string treeDesc = ss.str();
	return treeDesc;
}

// =========================================================================
//	Tree::writePaleoPLineage
// =========================================================================
void Tree::writePaleoPLineage(Node *p, stringstream &ss) {
	
	if (p != NULL){
		if(p->getLdes() == NULL) // tip
			ss << p->getTipName();
		else{
			ss << "(";
			writePaleoPLineage(p->getLdes(), ss);
			ss << makeFossPosString(p->getLdes());
			ss << ":" << p->getLdes()->getBranchLen();
			ss << ",";
			writePaleoPLineage(p->getRdes(), ss);
			ss << makeFossPosString(p->getRdes());
			ss << ":" << p->getRdes()->getBranchLen();
			ss << ")";
		}
	}
}

// =========================================================================
//	Tree::makeFossPosString
// =========================================================================
string Tree::makeFossPosString(Node *p) {
	
	stringstream ss;
	double fossDist = p->getNearestFAge();
	if(fossDist < 0.0)
		fossDist = totETL;
	int nF = p->getNumPPFossEvents();
	ss << "[&n_foss=" << nF;
	if(nF > 0){
		ss << ",fos_arr={" << p->getFossPosN(0);
		for(int i=1; i<nF; i++){
			ss << "," << p->getFossPosN(i);
		}
		ss << "}";

		ss << ",fos_tru_prob={" << p->getFossSampProbN(0);
		for(int i=1; i<nF; i++){
			ss << "," << p->getFossSampProbN(i);
		}
		ss << "}";
		
		Fossil *f = p->getFossilN(0);
		ss << ",fos_s_prob={" << f->getProbSamp();
		for(int i=1; i<nF; i++){
			f = p->getFossilN(i);
			ss << "," << f->getProbSamp();
		}
		ss << "}";
		
		int ns = p->getNumSampledMyFossils();
		if(ns > 0){
			double *sampFosArr = new double[ns];
			p->setSampFossAgesArray(sampFosArr);
			ss << ",fos_samp_bool={" << sampFosArr[0];
			for(int i=1; i<ns; i++)
				ss << "," << sampFosArr[i];
			ss << "}";
			delete [] sampFosArr;
		}
		
		
	}
	if(p->getIfNodeHasFossils()){
		vector<Fossil *> fF = p->getFoundFossilVector();
		ss << ",found_fos={" << fF[0]->getAge();
		for(int j=1; j<fF.size(); j++){
			double a = fF[j]->getAge();
			ss << "," << a;
		}
		ss << "}";
	}
	ss << ",anc_f_age=" << p->getAFAge() << ",ldes_f_age=" << p->getLDFAge() << ",rdes_f_age=" << p->getRDFAge();
	ss << ",near_dist=" << p->getNearFDist() << ",near_age=" << p->getNearestFAge();
	ss << ",f_dist=" << fossDist << "]";
	
	return ss.str();
}


// =========================================================================
//	Tree::evolveFossSamplingProbArr
// =========================================================================
void Tree::evolveFossSamplingProbArr(void) {
	
	int nf = 0;
	extantRoot->setSampPaleoProb(0.5);
	extantRoot->setIsExtantRoot(true);
	double rootD = extantRoot->getNodeDepth();
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			double ancFossPr = p->getAnc()->getSampPaleoProb();
			if(ancFossPr == 0.0){
				cout << "   Paleo prob == 0.0   " << endl;
			}
			double nBrLen = p->getBranchLen();
			int numEvents = p->getNumPPFossEvents();
			double newP = -1.0;
			if(numEvents > 0){
				nf += numEvents;
				p->initializeFossSampProb(numEvents);
				double brTime = 0.0;
				double oldP = ancFossPr;
				double oldTime = 0.0;
				for(int j=0; j<numEvents; j++){
					double fosPos = p->getFossPosN(j);
					brTime = fosPos - oldTime;
					newP = getNewFossSampProbability(oldP, brTime/rootD);
					p->setFossSampProbN(j, newP);
					oldP = newP;
					oldTime = fosPos;
				}
				brTime = nBrLen - oldTime;
				newP = getNewFossSampProbability(oldP, brTime/rootD);
				p->setSampPaleoProb(newP);
			}
			else{
				newP = getNewFossSampProbability(ancFossPr, nBrLen/rootD);
				p->setSampPaleoProb(newP);				
			}
			
		}
	}
	
	
	int counter = 0;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			p->initializeFossList(counter, includeRootFoss);
		}
	}
	numTotalFossilEvents = nf;
	if(includeRootFoss){
		setUpPaleoFossSamplingWithRoot();
	}
	else
		setUpPaleoFossSampling();
	
	
}

// =========================================================================
//	Tree::evolveFossRandomSamplingProbArr
// =========================================================================
void Tree::evolveFossRandomSamplingProbArr(void) {
	
	int nf = 0;
	extantRoot->setSampPaleoProb(0.5);
	extantRoot->setIsExtantRoot(true);
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			int numEvents = p->getNumPPFossEvents();
			if(numEvents > 0){
				nf += numEvents;
				p->initializeFossSampProb(numEvents);
				for(int j=0; j<numEvents; j++){
					p->setFossSampProbN(j, 1.0);
				}
				p->setSampPaleoProb(1.0);
			}
			else{
				p->setSampPaleoProb(1.0);				
			}
			
		}
	}
	
	int counter = 0;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			p->initializeFossList(counter, includeRootFoss);
		}
	}
	numTotalFossilEvents = nf;
	if(includeRootFoss){
		setUpPaleoFossSamplingWithRoot();
	}
	else
		setUpPaleoFossSampling();
}

// =========================================================================
//	Tree::getNewFossSampProbability
// =========================================================================
double Tree::getNewFossSampProbability(double exM, double bl) {
	double lnVar = kappaNuSamp * bl;
	double stDev = sqrt(lnVar);
	double mean = log(exM) - (lnVar * 0.5);//
	double newP = rando->logNormalRv(mean, stDev);
	return newP;
	
}


// =========================================================================
//	Tree::setUpPaleoFossSampling
// =========================================================================
void Tree::setUpPaleoFossSampling(void) {
	
	int sampleNum = paleoFossSampleNum;
	int nf = numTotalFossilEvents;
	if(nf < sampleNum)
		sampleNum = nf;
	double *fRates = new double[nf];
	Fossil *fList[nf];
	int counter = 0;
	double sumRate = 0.0;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot && p->getAnc() != extantRoot){
			int numEvents = p->getNumPPFossEvents();
			if(numEvents > 0){
				for(int j=0; j<numEvents; j++){
					fList[counter] = p->getFossilN(j);
					fRates[counter] = fList[counter]->getProbSamp();
					sumRate += fRates[counter];
					counter++;
				}
			}
		}
	}
	
	double *fProbs = new double[nf];
	
	Fossil *sampFoss[sampleNum];
	counter = 0;
	for(int f=0; f<sampleNum; f++){
		for(int i=0; i<nf; i++){
			double pr = fRates[i] / sumRate;
			if(i > 0)
				fProbs[i] = pr + fProbs[i-1];
			else
				fProbs[i] = pr;
		}
		if(fProbs[nf-1] < 1.0){
			for(int i=0; i<nf; i++)
				fProbs[i] /= fProbs[nf-1];
		}
		double ran = rando->uniformRv();
		int outcome = 0;
		while(ran > fProbs[outcome])
			outcome++;
		sampFoss[counter] = fList[outcome];
		counter++;
		sumRate -= fRates[outcome];
		fRates[outcome] = 0.0;
	}
	
	delete [] fRates;
	delete [] fProbs;
	
	vector<Node *> calNodes;
	
	for(int i=0; i<sampleNum; i++){
		Fossil *f = sampFoss[i];
		Node *p = f->getDecN()->getAnc();
		int flag = p->getFlag();
		while(flag != 2){
			p = p->getAnc();
			flag = p->getFlag();
		}
		p->addFossilToFoundVec(f);
		f->setIsFosSampled(true);
		calNodes.push_back(p);
	}
	
	
	setCalibratedFossilNodes();
	checkPaleoFossCompatibility();
	createCalibrationString();
	getAllFossilDistAncInfo();
	//cout << calibDataString;
}


// =========================================================================
//	Tree::setUpPaleoFossSamplingWithRoot
// =========================================================================
void Tree::setUpPaleoFossSamplingWithRoot(void) {
	
	int sampleNum = paleoFossSampleNum;
	int nf = numTotalFossilEvents;
	if(nf < sampleNum)
		sampleNum = nf;
	double *fRates = new double[nf];
	Fossil *fList[nf];
	int counter = 0;
	double sumRate = 0.0;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			int numEvents = p->getNumPPFossEvents();
			if(numEvents > 0){
				for(int j=0; j<numEvents; j++){
					fList[counter] = p->getFossilN(j);
					fRates[counter] = fList[counter]->getProbSamp();
					sumRate += fRates[counter];
					counter++;
				}
			}
		}
	}
	
	double *fProbs = new double[nf];
	
	Fossil *sampFoss[sampleNum];
	counter = 0;
	for(int f=0; f<sampleNum; f++){
		for(int i=0; i<nf; i++){
			double pr = fRates[i] / sumRate;
			if(i > 0)
				fProbs[i] = pr + fProbs[i-1];
			else
				fProbs[i] = pr;
		}
		if(fProbs[nf-1] < 1.0){
			for(int i=0; i<nf; i++)
				fProbs[i] /= fProbs[nf-1];
		}
		double ran = rando->uniformRv();
		int outcome = 0;
		while(ran > fProbs[outcome])
			outcome++;
		sampFoss[counter] = fList[outcome];
		counter++;
		sumRate -= fRates[outcome];
		fRates[outcome] = 0.0;
	}
	
	delete [] fRates;
	delete [] fProbs;
	
	vector<Node *> calNodes;
	
	for(int i=0; i<sampleNum; i++){
		Fossil *f = sampFoss[i];
		Node *p = f->getDecN()->getAnc();
		int flag = p->getFlag();
		while(flag != 2){
			p = p->getAnc();
			flag = p->getFlag();
		}
		p->addFossilToFoundVec(f);
		f->setIsFosSampled(true);
		calNodes.push_back(p);
	}
	
	setCalibratedFossilNodesWithRoot();
	checkPaleoFossCompatibilityWithRoot();
	createCalibrationStringFossRoot();
	getAllFossilDistAncInfo();
}

// =========================================================================
//	Tree::findLeftLeafOfNode
// =========================================================================
string Tree::findExtantLeftLeafOfNode(Node *p) {
	
	Node *d = p->getLdes();
	while(d->getIsTip() == false){
		if(!d->getLdes()->getIsExtant() && d->getLdes()->getIsTip())
			d = d->getRdes();
		else if(d->getLdes()->getFlag() < 1)
			d = d->getRdes();
		else
			d = d->getLdes();
			
	}
	return d->getTipName();
}	

// =========================================================================
//	Tree::findRightLeafOfNode
// =========================================================================
string Tree::findExtantRightLeafOfNode(Node *p) {
	
	Node *d = p->getRdes();
	while(d->getIsTip() == false){
		if(!d->getLdes()->getIsExtant() && d->getLdes()->getIsTip())
			d = d->getRdes();
		else if(d->getLdes()->getFlag() < 1)
			d = d->getRdes();
		else
			d = d->getLdes();
		
	}
	return d->getTipName();
}	


// =========================================================================
//	Tree::createCalibrationString
// =========================================================================
void Tree::createCalibrationString(void) {
	
	stringstream ss, ssFin, ssUnS;
	int nGoodFoss = 1;
	double tDepth = extantRoot->getNodeDepth();
	double minAge = extantRoot->getCalibFossilAge();
	double maxAge = tDepth + (tDepth * 0.1);
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p->getIfNodeHasFossils()){
			double myCal = p->getCalibFossilAge();
			
			
			string tL = findExtantLeftLeafOfNode(p);
			string tR = findExtantRightLeafOfNode(p);
			vector<Fossil *> fF = p->getFoundFossilVector();
			
			for(int j=0; j<fF.size(); j++){
				double fAge = fF[j]->getAge();
				if(fF[j]->getIsCalibrating() == true){
					nGoodFoss++;
					ss << "\n(" << tL << "," << tR << ")\t";
					if(fAge != myCal){
						cout << "fAge = " << fAge << " != myCal = " << myCal << endl;
					}
					ss << p->getNodeDepth() << "\t" << myCal << "\t-\t1";
				}
				else{
					ssUnS << "\n(" << tL << "," << tR << ")\t" << p->getNodeDepth() << "\t" << fAge << "\t-\t0";
				}
			}
		}
		
	}
	ssFin << numTotalFossilEvents << "\n";
	ssFin << "(root)\t" << tDepth << "\t" << minAge << "\t" << maxAge << "\t1";
	calibDataString = ssFin.str() + ss.str() + ssUnS.str();
	
}	

// =========================================================================
//	Tree::createCalibrationStringFossRoot
// =========================================================================
void Tree::createCalibrationStringFossRoot(void) {
	
	stringstream ss, ssFin, ssUnS;
	int nGoodFoss = 1;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p->getIfNodeHasFossils()){
			double myCal = p->getCalibFossilAge();
			
			
			string tL = findExtantLeftLeafOfNode(p);
			string tR = findExtantRightLeafOfNode(p);
			vector<Fossil *> fF = p->getFoundFossilVector();

			for(int j=0; j<fF.size(); j++){
				double fAge = fF[j]->getAge();
				double sProb = fF[j]->getProbSamp();
				if(fF[j]->getIsCalibrating() == true){
					nGoodFoss++;
					if(p == extantRoot)
						ss << "\n(root)\t";
					else
						ss << "\n(" << tL << "," << tR << ")\t";
					if(fAge != myCal){
						cout << "fAge = " << fAge << " != myCal = " << myCal << endl;
					}
					ss << p->getNodeDepth() << "\t" << myCal << "\t" << sProb << "\t1";
				}
				else{
					ssUnS << "\n(" << tL << "," << tR << ")\t" << p->getNodeDepth() << "\t" << fAge << "\t" << sProb << "\t0";
				}
			}
		}
		
	}
	ssFin << numTotalFossilEvents;
	calibDataString = ssFin.str() + ss.str() + ssUnS.str();
	
}	

// =========================================================================
//	Tree::setCalibratedFossilNodes
// =========================================================================
void Tree::setCalibratedFossilNodes(void) {
	
	int nGoodFoss = 1;
	double tDepth = extantRoot->getNodeDepth();
	double minAge = tDepth - (tDepth * 0.1);
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p->getIfNodeHasFossils()){
			nGoodFoss++;
			double oldest = 0.0;
			vector<Fossil *> fF = p->getFoundFossilVector();
			Fossil *calF;
			for(int j=0; j<fF.size(); j++){
				fF[j]->setIsCalibrating(false);
				double a = fF[j]->getAge();
				if(a > oldest){
					oldest = a;
					calF = fF[j];
				}
			}
			calF->setIsCalibrating(true);
			p->setCalibFossilAge(oldest);
			if(oldest > minAge)
				minAge = oldest;
		}
	}
	extantRoot->setCalibFossilAge(minAge);
	
	
}	


// =========================================================================
//	Tree::setCalibratedFossilNodesWithRoot
// =========================================================================
void Tree::setCalibratedFossilNodesWithRoot(void) {
	
	int nGoodFoss = 1;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p->getIfNodeHasFossils()){
			nGoodFoss++;
			double oldest = 0.0;
			vector<Fossil *> fF = p->getFoundFossilVector();
			Fossil *calF;
			for(int j=0; j<fF.size(); j++){
				fF[j]->setIsCalibrating(false);
				double a = fF[j]->getAge();
				if(a > oldest){
					oldest = a;
					calF = fF[j];
				}
			}
			calF->setIsCalibrating(true);
			p->setCalibFossilAge(oldest);
			
		}
	}
	
	
}	

// =========================================================================
//	Tree::checkPaleoFossCompatibility
// =========================================================================
int Tree::checkPaleoFossCompatibility(){
	
	int ncals = 0;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extDownPassSeq[i];
		if(p->getIfNodeHasFossils() && p != extantRoot){
			Node *q = p->getAnc();
			double myMinAge = p->getCalibFossilAge();
			double ancMinAge;
			do{
				if(q->getIfNodeHasFossils()){
					ancMinAge = q->getCalibFossilAge();
					if(myMinAge > ancMinAge){
						q->setCalibFossilAge(-1.0);
						q->falseFossilCalibrating();
					}
				}
				q = q->getAnc();
			}while(q != NULL);
		}
	}
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedNode())
			ncals++;
	}
	return ncals;
}	

// =========================================================================
//	Tree::checkPaleoFossCompatibilityWithRoot
// =========================================================================
int Tree::checkPaleoFossCompatibilityWithRoot(){
	
	int ncals = 0;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extDownPassSeq[i];
		if(p->getIfNodeHasFossils()){
			Node *q = p->getAnc();
			double myMinAge = p->getCalibFossilAge();
			double ancMinAge;
			if(q != NULL){
				do{
					if(q->getIfNodeHasFossils()){
						ancMinAge = q->getCalibFossilAge();
						if(myMinAge > ancMinAge){
							q->setCalibFossilAge(-1.0);
							q->falseFossilCalibrating();
						}
					}
					q = q->getAnc();
				}while(q != NULL);
			}
		}
	}
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedNode())
			ncals++;
	}
	return ncals;
}	

// =========================================================================
//	Tree::getAllFossilDistAncInfo
// =========================================================================
void Tree::getAllFossilDistAncInfo(void) {
	
	allFDString = "";
	stringstream ss;
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extUpPassSeq[i];
		if(p != extantRoot){
			int numEvents = p->getNumPPFossEvents();
			if(numEvents > 0){
				for(int j=0; j<numEvents; j++){
					Fossil *f = p->getFossilN(j);
					
					double ancND = f->getAncN()->getNodeDepth();
					double fosAge = f->getAge();
					
					double dist = ancND - fosAge;
					ss << dist << "\n";

				}
			}
		}
	}
	allFDString = ss.str();
	//cout << numTotalFossilEvents << endl;
}

// =========================================================================
//	Tree::calcNumAncDecFossilsExtantTree
// =========================================================================
void Tree::calcNumAncDecFossilsExtantTree(void) {
	
	//setExtantFlagsInFossilTree();
	exLinEvents = 0;
	int n = 0;
	
	int sumLins = 0;
	int extinctLs = 0;
	int living = 0;
	double tXTL = 0.0;
	
	for(int i=0; i<numExtPPTreeNodes; i++){
		Node *p = extDownPassSeq[i];
		int nFlag = p->getFlag();
		if(nFlag > 0){
			int nEvents = p->getNumPPFossEvents();
			exLinEvents += nEvents;
			n += nEvents;
			sumLins++;
			living++;
		}
		else{
			int nEvents = p->getNumPPFossEvents();
			n += nEvents;
			sumLins++;
			extinctLs++;
			tXTL += p->getBranchLen();
		}
	}
	numTotalFossilEvents = n;
	totalExtinctTL = tXTL;
	nExtinctLineages = extinctLs;
	extantRootDepth = extantRoot->getNodeDepth();
	
}

