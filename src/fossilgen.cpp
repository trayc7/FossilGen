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
#include <time.h>
#include <string.h>
#include "engine.h"

using namespace std;

void printHelp();
void printSettings(string of, int n, int r, double sc, double kapnu, double br, double dr, int ba, int da,
				   double bp, double dp);

void printHelp(){
	cout << "\tHere are the available options that you can change (default values are in []):\n";
	cout << "\t\t-out : output file name prefix **\n";
	cout << "\t\t-n   : number of extant taxa [= 100]\n";
	cout << "\t\t-r   : number of replicates [= 10]\n";
	cout << "\t\t-sc  : tree scale [= 1.0]\n";
	cout << "\t\t-br  : birth rate [= 0.5]\n";
	cout << "\t\t-dr  : death rate [= 0.2]\n";
	cout << "\t\t-s1  : seed 1 (use this if you only pass in one seed) \n";
	cout << "\t\t-s1  : seed 2 \n";
	cout << "\t\t-pp  : paleo simulation psi rate [= 0.1] \n";
	cout << "\t\t-sga : scale tree from gamma with alpha [= 4.0] \n";
	cout << "\t\t-sgr : scale tree from gamma with rate [= 8.0] \n";
	cout << "\t\t-fsp : number of fossils to sample\n";
	cout << "\t\t-rnd : randomly sample all fossils\n";
	cout << "\t\t-kap : sample fossils in proportion to sampling rate with variance parameter [= 0.005]\n";

}

void printSettings(string of, int n, int r, double sc, double kapnu, double br, double dr, int ba, int da,
				   double bp, double dp, double psi){
	

	cout << "\t\toutput file name prefix    = " << of << "\n";
	cout << "\t\tnumber of extant taxa      = " << n << "\n";
	cout << "\t\tnumber of reps             = " << r << "\n";
	cout << "\t\tbirth rate                 = " << br << "\n";
	cout << "\t\tdeath rate                 = " << dr << "\n";
	cout << "\t\tfossilization rate         = " << psi << "\n";
	cout << "\t\tpreservation ac parameter  = " << kapnu << "\n";
}

int main(int p, char *argStr[])
{
	PhyEngine *runPhySim;
	string setFileName = "";
	if(p == 1){
		printHelp();
		return 0;
	}
	else {
		string outFN = "";
		int modelType=3, birthAlph=100, deathAlph=100, reps=10, ntaxa=100, seed1=-1, seed2=-1, fossSamp=-1;
		double birthR=0.5, deathR=0.2, birthPrior=3.0, deathPrior=3.0, treeScale=1.0,
				phiRate=0.1, scaleAlph=4.0, scaleRate=8.0, kappaNu = 0.005;
		bool palP = true, doGammScale=false, includeRoot=true, randomFossSample=true, sampAllFoss=false;
		for (int i = 1; i < p; i++){
			char *curArg = argStr[i];
			if(strlen(curArg) > 1 && curArg[0] == '-'){
				if(!strcmp(curArg, "-out"))
					outFN = argStr[i+1];
				else if(!strcmp(curArg, "-m"))
					modelType = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-br"))
					birthR = atof(argStr[i+1]);
				else if(!strcmp(curArg, "-dr"))
					deathR = atof(argStr[i+1]);
				else if(!strcmp(curArg, "-n"))
					ntaxa = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-ba"))
					birthAlph = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-da"))
					deathAlph = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-bp"))
					birthPrior = atof(argStr[i+1]);
				else if(!strcmp(curArg, "-dp"))
					deathPrior = atof(argStr[i+1]);
				else if(!strcmp(curArg, "-sc"))
					treeScale = atof(argStr[i+1]);
				else if(!strcmp(curArg, "-r"))
					reps = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-set"))
					setFileName = argStr[i+1];
				else if(!strcmp(curArg, "-h")){
					printHelp();
					return 0;
				}
				else if(!strcmp(curArg, "-pp")){
					palP = true;
					phiRate = atof(argStr[i+1]);
				}
				else if(!strcmp(curArg, "-kap")){
					kappaNu = atof(argStr[i+1]);
					randomFossSample = false;
				}
				else if(!strcmp(curArg, "-sga")){
					scaleAlph = atof(argStr[i+1]);
					doGammScale=true;
				}
				else if(!strcmp(curArg, "-sgr"))
					scaleRate = atof(argStr[i+1]);
				else if(!strcmp(curArg, "-rnd"))
					randomFossSample = true;
				else if(!strcmp(curArg, "-fsp"))
					fossSamp = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-allf")){ // Need this to plot all
					sampAllFoss = true;
					fossSamp = 2;
				}
				else if(!strcmp(curArg, "-s1"))
					seed1 = atoi(argStr[i+1]);
				else if(!strcmp(curArg, "-s2"))
					seed2 = atoi(argStr[i+1]);
				else {
					cout << "\n############################ !!! ###########################\n";
					cout << "\n\n\tPerhaps you mis-typed something, here are the \n\tavailable options:\n";
					printHelp();
					cout << "\n############################ !!! ###########################\n";
					return 0;
				}
			}
		}
		if(setFileName.empty()){
			printSettings(outFN, ntaxa, reps, treeScale, kappaNu, birthR, deathR, birthAlph, deathAlph,
						  birthPrior, deathPrior, phiRate);
			
			runPhySim = new PhyEngine(outFN, modelType, birthR, deathR, ntaxa, birthAlph, deathAlph, birthPrior, 
									  deathPrior, treeScale, reps, seed1, seed2, palP, phiRate, scaleAlph, 
									  scaleRate, doGammScale, fossSamp, includeRoot, randomFossSample, kappaNu,
									  sampAllFoss);
		}
	}
	
	runPhySim->doRunRun();
	
	delete runPhySim;
	return 0;
}
