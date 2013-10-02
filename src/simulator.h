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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <string>
#include "MbRandom.h"


class Tree;

/*----------------------------------------------
 class Simulator
 ------------------------------------------------*/

class Simulator{
	private:
		Tree							*simt;
		Tree							*extantTree;
		Tree							*fossilTree;
		int								modeltype;
		int								extantStop, gsaExtStop;
		double							treescale;
		double							stBirthRt, stDeathRt;
		int								birthSampAlpha, deathSampAlpha; 
		double							bPriorAlpha, dPriorAlpha; 
		MbRandom						*rando;
		double							currentSimTime;
		std::vector<Tree*>				gsaTreess;
		bool							includeRootFossils, allFossilsSamp;

		
	public:
						Simulator(MbRandom *rndp, int mt, int nt, double sc, double br, double dr, 
								  int ba, int da, double bp, double dp, bool inRt, bool allF);
						~Simulator();
		bool			simulateTree();
		void			initializeSimulator();
		bool			bdSim();
		bool			varSim();
		bool			gsaBDSim();
		bool			gsaVarSim();
		bool			czechStop();
		void			processSimulation(bool scT);
		void			processGSASimulation();
		std::string		getWholeTreeStringSim();
		std::string		getExtantTreeStringSim();
		std::string		getRawCalScaledTreeStringSim();
		double			getExtantTreeDepth();
		void			setExtGSAStop(int extStop, int gsaStop);


		void			doPaleoSimStuff(std::ofstream &dOut, double phiF, int numFossToSamp, bool rndFosSamp, double kpnu);
		std::string		getPaleoSimNewickTreeDesc(void);
		std::string		getPaleoSimCalibrationDataString(void);
		std::string		getAllPaleoFossDistDataString(void);
		std::string		getTreeNumFossilExtantLinFossilString(void);
		
};	


#endif