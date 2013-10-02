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

#ifndef TREE_H
#define TREE_H

#include <set>
#include <string>
#include <vector>
#include "MbRandom.h"

class Node;
/*----------------------------------------------
 class Fossil
 ------------------------------------------------*/
class Fossil{
	
	private:
		Node			*decN;
		Node			*ancN;
		double			age;
		double			probSamp;
		int				index;
		double			trueProb;
		bool			isSampled;
		bool			isCalib;
		
	public:
						Fossil(void);
		double			getAge(void) { return age; }
		void			setAge(double a) { age = a; }
		double			getProbSamp(void) { return probSamp; }
		void			setProbSamp(double s) { probSamp = s; }
		Node*			getDecN(void) { return decN; }
		void			setDecN(Node *p) { decN = p; }
		Node*			getAncN(void) { return ancN; }
		void			setAncN(Node *p) { ancN = p; }
		int				getFossIndex(void) { return index; }
		void			setFossIndex(int i) { index = i; }
		double			getTrueProbSamp(void) { return trueProb; }
		void			setTrueProbSamp(double s) { trueProb = s; }
		bool			getIsFosSampled(void) { return isSampled; }
		void			setIsFosSampled(bool b) { isSampled = b; }
		bool			getIsCalibrating(void) { return isCalib; }
		void			setIsCalibrating(bool b) { isCalib = b; }

};


/*----------------------------------------------
 class Node
 ------------------------------------------------*/
class Node{
	private:
		Node			*ldes;
		Node			*rdes;
		Node			*anc;
		Node			*sib;
		int				indx;
		std::string		tname;
		bool			isTip;
		double			branchLen;
		double			nodeDepth, birthTime, deathTime;
		double			birthRate, deathRate;
		int				flag;
		bool			isRoot, isExtant;
		double			prSampFoss;
		bool			isIGFossil, isFoundFossil, isCalibratedNode;
		int				numDecTaxa;
		double			calibMinAge;
		double			*randTestCals;
		

		double			*fossPositions;
		double			*fossSampProbs;
		int				numFossEvents;
		int				numAllDecTaxa;
		double			lDFD, rDFD, aFD, nearFoss;
		double			lDFTime, rDFTime, aFTime, nearTime;
		Node			*lDirF, *rDirF, *aDirF;
		double			sampPaleoProb;
		Fossil			*myFossils;
		std::vector<Fossil *>	foundFossils;
		bool			isExtRoot;
		double			calFossAge;
		
	public:
						Node(void);
						~Node(void);
		void			nullNode(void);
		void			nullPPFD(void);
		void			initializeRandTestCals(int n) { randTestCals = new double[n]; }
		double			*getRandTestCalsAll() { return randTestCals; }
		double			getRandTestCalN(int n) { return randTestCals[n]; }
		void			setRandTestCalN(int n, double v) { randTestCals[n] = v; }
		Node*			getLdes(void) { return ldes; }
		Node*			getSib(void) { return sib; }
		Node*			getRdes(void) { return rdes; }
		Node*			getAnc(void) { return anc; }
		int				getIndex(void) { return indx; }
		double			getBranchLen(void) { return branchLen; }
		std::string		getTipName(void) { return tname; }
		double			getNBirthRate(void) { return birthRate; }
		double			getNDeathRate(void) { return deathRate; }
		int				getFlag(void) { return flag; }
		bool			getIsTip() { return isTip; }
		bool			getIsExtant() { return isExtant; }
		bool			getIsRoot() { return isRoot; }
		double			getNodeDepth() { return nodeDepth; }
		double			getNodeBirthTime() { return birthTime; }
		double			getNodeDeathTime() { return deathTime; }
		double			getPrSampFossil() { return prSampFoss;}
		bool			getIsIGFossil() { return isIGFossil; }
		bool			getIsFoundFossil() { return isFoundFossil; }
		bool			getIsCalibratedNode() { return isCalibratedNode; }
		int				getNumDecExtantTax() { return numDecTaxa; }
		double			getCalibrationMinAge() { return calibMinAge; }
		void			setIndex(int i) { indx = i; }
		void			setLdes(Node *p) { ldes = p; }
		void			setSib(Node *p) { sib = p; }
		void			setRdes(Node *p) { rdes = p; }
		void			setAnc(Node *p) { anc = p; }
		void			setTipName(std::string n) { tname = n; }
		void			setIsTip(bool b) { isTip = b; }
		void			setBranchLen(double v) { branchLen = v; }
		void			setNBirthRate(double v) { birthRate = v; }
		void			setNDeathRate(double v) { deathRate = v; }
		void			setFlag(int f) { flag = f; }
		void			setIsRoot(bool b) { isRoot = b; }
		void			setNodeDepth(double v) { nodeDepth = v; }
		void			setNodeBirthTime(double v) { birthTime = v; }
		void			setNodeDeathTime(double v) { deathTime = v; }
		void			setIsExtant(bool b) { isExtant = b; }
		void			setPrSampFossil(double p) { prSampFoss = p;}
		void			setIsIGFossil(bool b) { isIGFossil = b; }
		void			setIsFoundFossil(bool b) { isFoundFossil = b; }
		void			setIsCalibratedNode(bool b) { isCalibratedNode = b; }
		void			setNumDecExtantTax(int i) { numDecTaxa = i; }
		void			setCalibrationMinAge(double v) { calibMinAge = v; }
	
		void			initializeFossPos(int n) { fossPositions = new double[n]; }
		void			initializeFossSampProb(int n) { fossSampProbs = new double[n]; }
		double			*getFossPosAll() { return fossPositions; }
		double			getFossPosN(int n) { return fossPositions[n]; }
		void			setFossPosN(int n, double v) { fossPositions[n] = v; }
		double			getFossSampProbN(int n) { return fossSampProbs[n]; }
		void			setFossSampProbN(int n, double v) { fossSampProbs[n] = v; }
		int				getNumPPFossEvents() { return numFossEvents; }
		void			setNumPPFossEvents(int n) { numFossEvents = n; }
		int				getNumAllDecExtantTax() { return numAllDecTaxa; }
		void			setNumAllDecExtantTax(int i) { numAllDecTaxa = i; }
		double			getLDFDist() { return lDFD; }
		void			setLDFDist(double d) { lDFD = d; }
		double			getRDFDist() { return rDFD; }
		void			setRDFDist(double d) { rDFD = d; }
		double			getAFDist() { return aFD; }
		void			setAFDist(double d) { aFD = d; }
		double			getNearFDist() { return nearFoss; }
		void			setNearFDist(double d) { nearFoss = d; }
		Node*			getLDirFoss(void) { return lDirF; }
		Node*			getRDirFoss(void) { return rDirF; }
		Node*			getADirFoss(void) { return aDirF; }
		void			setLDirFoss(Node *n) { lDirF = n; }
		void			setRDirFoss(Node *n) { rDirF = n; }
		void			setADirFoss(Node *n) { aDirF = n; }
		double			getAgeOfFossilAtPosN(int n) { return nodeDepth + branchLen - fossPositions[n]; }
		double			getLDFAge() { return lDFTime; }
		void			setLDFAge(double d) { lDFTime = d; }
		double			getRDFAge() { return rDFTime; }
		void			setRDFAge(double d) { rDFTime = d; }
		double			getAFAge() { return aFTime; }
		void			setAFAge(double d) { aFTime = d; }
		double			getNearestFAge() { return nearTime; }
		void			setNearestFAge(double d) { nearTime = d; }
		double			getSampPaleoProb() { return sampPaleoProb; }
		void			setSampPaleoProb(double d) { sampPaleoProb = d; }
		void			initializeFossList(int &ix, bool inclRoot);
		double			getAgeForFossIndex(int i);
		Fossil			*getFossilN(int n) { return &myFossils[n]; }
		void			addFossilToFoundVec(Fossil *f) { foundFossils.push_back(f); }
		bool			getIfNodeHasFossils(void) { return !foundFossils.empty(); }
		std::vector<Fossil *>	getFoundFossilVector(void) { return foundFossils; }
		void			setIsExtantRoot(bool b) { isExtRoot = b; }
		bool			getIsExtantRoot(void) { return isExtRoot; }
		int				getNumSampledMyFossils(void);
		void			setSampFossAgesArray(double *a);
		double			getCalibFossilAge(void) { return calFossAge; }
		void			setCalibFossilAge(double d) { calFossAge = d; }
		void			falseFossilCalibrating(void);



};

/*----------------------------------------------
 class Tree
 ------------------------------------------------*/

class Tree{
	
	private:
		Node			*root;
		int				ntax, numNodes, numTotalTips;
		int				treeindex;
		Node			*nodes;
		
		Node			**upPassSeq;
		MbRandom		*rando;
		double			bRate0, dRate0;
		
		std::vector<Node *>	simulating;
		std::vector<Node *>	extants;
		std::vector<Node *>	inGFossils;
		
		int				allocated, available, capacity;
		int				numExtant, numExtinct;
		int				modelsim;
		
		int				birthSampAlpha, deathSampAlpha; 
		double			bPriorAlpha, dPriorAlpha; 
		
		double			currentTime;
		double			rawTreeDepth, scaleToVal;
		
		int				nSampT;
		
		Node			*extantRoot;
		Node			**extUpPassSeq;
		Node			**extDownPassSeq;
		int				numExtPPTreeNodes;
		int				numExtPPTreeTips;
		double			totETL;
		
		std::string		calibDataString;
		int				numTotalFossilEvents;
		int				paleoFossSampleNum;
		std::string		allFDString;
		bool			includeRootFoss, rndFossSamp;
		int				exLinEvents, nExtinctLineages;
		double			totalExtinctTL, extantRootDepth;
		double			kappaNuSamp;
		bool			sampAllFossils;
		

	
	public:
						Tree(MbRandom *p, double b0, double d0, int ext, double ct, int ba, int da,
							 double bp, double dp, bool inRt, bool allF);
						Tree(MbRandom *p, int nt);
						~Tree();
		
		
		void			initializeTree(double curTime);
		int				getNumExtantLineages() { return numExtant; }
		int				getNumExtinctLineages() { return numExtinct; }
		void			setNewLineageInfo(Node *par, Node *offsp, Node *sib);
		double			getTimeToNextERMEvent();
		double			getTimeToNextVarEvent();
		void			ermEvent(double ct);
		void			varEvent(double ct);
		void			setExtantPresentTime(double ct);
		void			lineageDeathEvent(Node *p, int idx);
		void			lineageBirthEvenet(Node *p, int idx);
		void			mutateLineage(Node *par, Node *dec);
		void			gammaMutateLineage(Node *offsp, double pb, double pd);
		void			zeroAllFlags();
		void			setWholeTreeFlags();
		void			setExtantTreeFlags();
		void			setGSATipTreeFlags();
		void			popExtantVec();
		void			setSampleFromFlags();
		void			setWholeTreeExtantRoot();
		void			setTreeBranchLengths();
		void			setTreeNodeDepths();
		void			setTreeTipNames();
		void			setRecTreeTipNames(Node *inN, int &extIt, int &tipIt);
		std::string		getWholeTreeString();
		void			reconstructTreeFromSimulation(Node *oRoot);
		int				reconstructLineageFromSimulation(Node *myN, Node *oN, int &nextIntN, int &tipcounter);
		void			passUp(Node *p, int *x);
		void			getRevDownPassSequence();
		std::string		getNewickTreeDesc();
		std::string		getIncWholeTreeNewick();
		void			writeTreeLineage(Node *p, std::stringstream &ss);
		Node			*getTreeRoot() { return root; }
		void			prepTreeForReconstruction();
		void			prepGSATreeForReconstruction();
		void			scaleTreeDepthToValue(double scVal);
		void			setModelSimType(int t) { modelsim = t; }
		std::string		getCalScaleTreeDesc();
		void			writeCalScaleTreeLineage(Node *p, std::stringstream &ss);
		std::string		findLeftLeafOfNode(Node *p);
		std::string		findRightLeafOfNode(Node *p);
		double			getTotalTreeDepth() { return rawTreeDepth; }
		
		Node*			getExtantTreeRoot(Node *p);
		std::string		getFossilNewickTreeDesc(void);
		void			writeFossilTreeLineage(Node *p, std::stringstream &ss);
		int				setFossilExtantTreeFlags();
		int				getTreeNumDecTaxa();
		int				getNodeNumDecTaxa(Node *p);
		void			setAllTreeNodeIDs();
		void			setTreeNodeIDs(Node *p, int &tipID, int &nID);
		void			setExtantFlagsInFossilTree();
		std::string		findLeftLeafOfFossilNode(Node *p);
		std::string		findRightLeafOfFossilNode(Node *p);
		int				getNumberOfCalibNodes();
		double			getNodeTimeLowerBound(Node *p);
		void			quickSort(double numbers[], int arrSize);
		void			qSort(double numbers[], int left, int right);
		void			setNSampT(int i) { nSampT = i; }
		void			setExtTreeExtRoot() { extantRoot = root; }
		void			setAllFossilSample(bool t) { sampAllFossils = t; }
		
		
		void			generatePaleoPriorFoss(std::ofstream &dOut, double phi, int sampNum, bool rfs, double kpnu);
		Node*			getNodeNearestPPFoss(Node *p, Node *a);
		void			getRevDPSeqExtRootTree(void);
		void			passExtTUp(Node *p, int *x);
		void			getDownPSeqExtRootTree(void);
		void			passExtTDown(Node *p, int *x);
		int				countNTipsInExtFosTree(void);
		int				setNumTipDec(Node *p);
		double			getTotalTreeLength(void);
		double			getNodesDistToClosestFoss(Node *p);
		void			nullNodePPFossInfo();
		void			setRecNodeFosPositions(void);
		double			getDistBWTwoNodes(Node *p, Node *q, int updown);
		std::string		getPaleoPNewickDesc(void);
		void			writePaleoPLineage(Node *p, std::stringstream &ss);
		std::string		makeFossPosString(Node *p);
		void			evolveFossSamplingProbArr(void);
		void			evolveFossRandomSamplingProbArr(void);
		double			getNewFossSampProbability(double exM, double bl);
		void			setUpPaleoFossSampling(void);
		void			setUpPaleoFossSamplingWithRoot(void);
		std::string		findExtantLeftLeafOfNode(Node *p);
		std::string		findExtantRightLeafOfNode(Node *p);
		void			createCalibrationString(void);
		void			createCalibrationStringFossRoot(void);
		std::string		getCalibrationDataString(void) { return calibDataString; }
		void			setCalibratedFossilNodes(void);
		void			setCalibratedFossilNodesWithRoot(void);
		int				checkPaleoFossCompatibility(void);
		int				checkPaleoFossCompatibilityWithRoot(void);
		void			getAllFossilDistAncInfo(void);
		std::string		getFossDistAncString(void){ return allFDString; }
		void			setIncludeRootFossSampBool(bool b) { includeRootFoss = b; }
		void			calcNumAncDecFossilsExtantTree(void);
		int				getNumExLineageFossils(void) { return exLinEvents; }
		int				getNumTotalFossilEvents(void) { return numTotalFossilEvents; }
		int				getNumExtinctLineagesOnTree(void) { return nExtinctLineages; }
		double			getTotalTreeLengthWholeT(void) { return totETL; }
		double			getTotalLenghtOfExtinctLineages(void) { return totalExtinctTL; }
		int				getNumExtPPTreeNodes(void) { return numExtPPTreeNodes; }
		double			getExtantRootDepth(void) { return extantRootDepth; }
};	







#endif