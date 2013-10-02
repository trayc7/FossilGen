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

#ifndef ENGINE_H
#define ENGINE_H

#include <vector>
#include "MbRandom.h"
 

/*----------------------------------------------
 class TreeInfo
 ------------------------------------------------*/
class TreeInfo{
	private:
		int							tInd;
		std::string					wholeT, extantTree, timeScaleTree;
		std::string					wholeFig, extantFig, fossilFig;
		std::vector<std::string>	gsaTree;
		double						treelen, treeness, avetiplen, treeDepth;
		
		std::string					pfCalibData, pfPaleoTreeFig;

		
	public:
							TreeInfo(int idx) : tInd(idx), treelen(0.0), treeness(0.0),
												avetiplen(0.0) {}
		std::string				getWholeTreeStringInfo() { return wholeT; }
		std::string				getExtantTreeStringInfo() { return extantTree; }
		std::string				getWholeFigStringInfo() { return wholeFig; }
		std::string				getExtantFigStringInfo() { return extantFig; }
		std::string				getTimeTreeStringInfo() { return timeScaleTree; }
		std::string				getFossilTreeStringInfo() { return fossilFig; }
		std::vector<std::string>				getGSATreeStringInfo() { return gsaTree; }
		double				getTreeLenInfo() { return treelen; }
		double				getTreenessInfo() { return treeness; }
		double				getAveTipLenInfo() { return avetiplen; }
		int					getTreeIDInfo() { return tInd; }
		double				getTreeDepthInfo() { return treeDepth; }
		
		void				setWholeTreeStringInfo(std::string ts) { wholeT = ts; }
		void				setExtantTreeStringInfo(std::string ts) { extantTree = ts; }
		void				setWholeFigStringInfo(std::string ts) { wholeFig = ts; }
		void				setExtantFigtringInfo(std::string ts) { extantFig = ts; }
		void				setTimeTreeStringInfo(std::string ts) { timeScaleTree = ts; }
		void				setFossilTreeStringInfo(std::string ts) { fossilFig = ts; }
		void				setGSATreeStringInfo(std::vector<std::string> ts) { gsaTree = ts; }
		void				setTreeLenInfo(double v) { treelen = v; }
		void				setTreenessInfo(double v) { treeness = v; }
		void				setAveTipLenInfo(double v) { avetiplen = v; }
		void				setTreeDepthInfo(double v) { treeDepth = v; }
		
		void				writeWholeTreeFileInfo(std::string ofp);
		void				writeExtantTreeFileInfo(std::string ofp);
		void				writeScaledDepthTreeFileInfo(std::string ofp);
		void				writeFossilTreeFileInfo(std::string ofp);
		void				writeGSATreeFileInfo(std::string ofp);
		
		std::string			getPaleoFossCalibData() { return pfCalibData; }
		void				setPaleoFossCalibData(std::string ts) { pfCalibData = ts; }
		std::string			getPaleoFossFigTreeString() { return pfPaleoTreeFig; }
		void				setPaleoFossFigTreeString(std::string ts) { pfPaleoTreeFig = ts; }
		void				writePaleoFossilTreeFileInfo(std::string ofp);


		
};


/*----------------------------------------------
 class PhyEngine
 ------------------------------------------------*/
class PhyEngine{
	private:
	
		std::string				outfilename;
		std::vector<TreeInfo*>	simtrees;
		
		int					modelType;
		int					seedset;
		int					extantStop, numTrees;
		double				treescale;
		bool				doScaleTree;
		
		MbRandom			rando;
		
		double				stBirthRt, stDeathRt;
		int					birthSampAlpha, deathSampAlpha; 
		double				bPriorAlpha, dPriorAlpha; 
		bool				doPaleoPrior, allFossils;
		double				fossRatePhi;
		double				scaleGammaAlpha, scaleGammaRate;
		bool				doScaleRanGamma, includeRootFossilSamp, rndFossSamp;
		int					numFossilSamples;
		double				kappaNu; 
		
		
	
	public:
							PhyEngine(std::string of, int mt, double br, double dr, 
									  int nt, int ba, int da, double bp, double dp, 
									  double sc, double rp, int sd1, 
									  int sd2, bool ppp, double frt, double sgal,
									  double sgrt, bool doSG, int fsamp, bool inRt,
									  bool rnsmp, double kpanu, bool allFos);
							~PhyEngine();
		void				doRunRun();
		void				doSimulatePaleoPriorFoss();
		void				writeTreeFiles();
		TreeInfo			*findTreeByID(int i);
		void				calcAverageRootAges();
		
};

#endif