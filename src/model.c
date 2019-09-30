/*
 *  MrBayes 3.1.2
 *
 *  copyright 2002-2005
 *
 *  John P. Huelsenbeck
 *  Section of Ecology, Behavior and Evolution
 *  Division of Biological Sciences
 *  University of California, San Diego
 *  La Jolla, CA 92093-0116
 *
 *  johnh@biomail.ucsd.edu
 *
 *	Fredrik Ronquist
 *  School of Computational Science
 *  Florida State University
 *  Tallahassee, FL 32306-4120
 *
 *  ronquist@csit.fsu.edu
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "mb.h"
#include "globals.h"
#include "bayes.h"
#include "model.h"
#include "command.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif


/* local prototypes */
int IsModelSame (int whichParam, int part1, int part2, int *isApplic1, int *isApplic2);
int NumActiveParts (void);


/* globals */
Model			modelParams[MAX_NUM_DIVS];              /* holds model params for up to MAX_NUM_DIVS partitions    */
int				numCurrentDivisions;          			/* number of partitions of data                            */
int				activeParams[NUM_LINKED][MAX_NUM_DIVS]; /* a table holding the parameter link status               */

/* local */
static int		numVars[MAX_NUM_DIVS], activeParts[MAX_NUM_DIVS], fromI, toJ, foundDash, foundComma, foundEqual, foundBeta, foundAaSetting, foundExp, modelIsFixed,
				linkNum, linkTable[NUM_LINKED][MAX_NUM_DIVS], tempLinkUnlink[NUM_LINKED][MAX_NUM_DIVS], tempLinkUnlinkVec[MAX_NUM_DIVS], tempNumStates, isNegative;
static MrBFlt	tempStateFreqs[200], tempAaModelPrs[10], tempNum[MAX_NUM_DIVS];
static char		colonPr[100];





int AreDoublesEqual (MrBFlt x, MrBFlt y, MrBFlt tol)

{

	if ((x - y) < -tol || (x - y) > tol)
		return (NO);
	else
		return (YES);
	
}





int CheckModel (void)

{
	
	int			a, b, d, i, j, k, m, n, ns, paramCount, nOfParam, isFirst, shouldPrint, isSame;
	char		tempName[100];

	/* set up table of linkages */
	if (SetUpModel () == ERROR)
		{
		MrBayesPrint ("%s   Problem setting up default links\n", spacer);
		return (ERROR);
		}

	MrBayesPrint ("%s   Model settings:\n\n", spacer);
	for (i=0; i<numCurrentDivisions; i++)
		{
		ns = 0;

		if (numCurrentDivisions > 1)
			MrBayesPrint ("%s      Settings for partition %d --\n", spacer, i+1);
		
		if (modelParams[i].dataType == DNA)
			{
			MrBayesPrint ("%s         Datatype  = DNA\n", spacer);
			ns = 4;
			}
		else if (modelParams[i].dataType == RNA)
			{
			MrBayesPrint ("%s         Datatype  = RNA\n", spacer);
			ns = 4;
			}
		else if (modelParams[i].dataType == PROTEIN)
			{
			MrBayesPrint ("%s         Datatype  = Protein\n", spacer);
			ns = 20;
			}
		else if (modelParams[i].dataType == RESTRICTION)
			{
			MrBayesPrint ("%s         Datatype  = Restriction\n", spacer);
			ns = 2;
			}
		else if (modelParams[i].dataType == STANDARD)
			{
			MrBayesPrint ("%s         Datatype  = Standard\n", spacer);
			ns = 10;
			}
		else if (modelParams[i].dataType == CONTINUOUS)
			{
			MrBayesPrint ("%s         Datatype  = Continuous\n", spacer);
			}
			
		if (modelParams[i].dataType == CONTINUOUS)
			{
			/* begin description of continuous models */
			  if (!strcmp(modelParams[i].brownCorPr, "Fixed") && AreDoublesEqual(modelParams[i].brownCorrFix, 0.0, ETA)==YES)
				MrBayesPrint ("%s         Model     = Independent Brownian motion\n", spacer);
			else
				MrBayesPrint ("%s         Model     = Correlated Brownian motion\n", spacer);
			/* end description of continuous models */
			}
		else
			{
			/* begin description of discrete models */
			if (!strcmp(modelParams[i].parsModel, "Yes"))
				{
				MrBayesPrint ("%s         Parsmodel = %s\n", spacer, modelParams[i].parsModel);
				}
			else
				{
				/* dna characters in this partition */
				if (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA)
					{
					/* general form of the rate matrix */ 
					MrBayesPrint ("%s         Nucmodel  = %s\n", spacer, modelParams[i].nucModel);
				
					/* constraints on rates of substitution */
					MrBayesPrint ("%s         Nst       = %s\n", spacer, modelParams[i].nst);
					if (!strcmp(modelParams[i].nst, "2"))
						{
						if (!strcmp(modelParams[i].tRatioPr,"Beta"))
							{
							MrBayesPrint ("%s                     Transition and transversion  rates, expressed\n", spacer);
							MrBayesPrint ("%s                     as proportions of the rate sum, have a\n", spacer);
							MrBayesPrint ("%s                     Beta(%1.2lf,%1.2lf) prior\n", spacer, modelParams[i].tRatioDir[0], modelParams[i].tRatioDir[1]);
							}
						else
							{
							MrBayesPrint ("%s                     Transition/transversion rate ratio is fixed to %1.2lf.\n", spacer, modelParams[i].tRatioFix);
							}
						}
					else if (!strcmp(modelParams[i].nst, "6"))
						{
						if (!strcmp(modelParams[i].revMatPr,"Dirichlet"))
							{
							MrBayesPrint ("%s                     Substitution rates, expressed as proportions\n", spacer);
							MrBayesPrint ("%s                     of the rate sum, have a Dirichlet prior\n", spacer);
							MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
								modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
								modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
							}
						else
							{
							MrBayesPrint ("%s                     Substitution rates are fixed to be \n", spacer);
							MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf).\n", spacer, modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2], modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);
							}
						}
					else if (!strcmp(modelParams[i].nst, "Mixed"))
						{
						if (!strcmp(modelParams[i].revMatPr,"Dirichlet"))
							{
							MrBayesPrint ("%s                     Substitution rates, expressed as proportions\n", spacer);
							MrBayesPrint ("%s                     of the rate sum, have a Dirichlet prior\n", spacer);
							MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
								modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
								modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
							}
						else
							{
							MrBayesPrint ("%s                     Substitution rates are fixed to be \n", spacer);
							MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf).\n", spacer, modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2], modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);
							}
						}
					
					if (!strcmp(modelParams[i].nucModel,"Codon"))
						{
						/* what is the distribution on the nonsyn./syn. rate ratio */
						if (!strcmp(modelParams[i].omegaVar, "Equal"))
							{
							if (!strcmp(modelParams[i].omegaPr,"Dirichlet"))
								{
								MrBayesPrint ("%s                     Nonsynonymous and synonymous rates, expressed\n", spacer);
								MrBayesPrint ("%s                     as proportions of the rate sum, have a\n", spacer);
								MrBayesPrint ("%s                     Dirichlet(%1.2lf,%1.2lf) prior\n", spacer, modelParams[i].omegaDir[0], modelParams[i].omegaDir[1]);
								}
							else
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio is fixed to %1.2lf.\n", spacer, modelParams[i].omegaFix);
								}
							}
						else if (!strcmp(modelParams[i].omegaVar, "Ny98"))
							{
							if (!strcmp(modelParams[i].ny98omega1pr, "Beta"))
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for purifying selection\n", spacer);
								MrBayesPrint ("%s                     (class 1) has a Beta(%1.2lf,%1.2lf) on the interval (0,1).\n", spacer, modelParams[i].ny98omega1Beta[0], modelParams[i].ny98omega1Beta[1]);
								}
							else
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for purifying selection\n", spacer);
								MrBayesPrint ("%s                     (class 1) is fixed to %1.2lf.\n", spacer, modelParams[i].ny98omega1Fixed);
								}
							MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for neutral selection\n", spacer);
							MrBayesPrint ("%s                     (class 2) is fixed to 1.0.\n", spacer);
							if (!strcmp(modelParams[i].ny98omega3pr, "Uniform"))
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for positive selection\n", spacer);
								MrBayesPrint ("%s                     is uniformly distributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].ny98omega3Uni[0], modelParams[i].ny98omega3Uni[1]);
								}
							else if (!strcmp(modelParams[i].ny98omega3pr, "Exponential"))
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for positive selection\n", spacer);
								MrBayesPrint ("%s                     is exponentially distributed with parameter (%1.2lf).\n", spacer, modelParams[i].ny98omega3Exp);
								}
							else
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for positive \n", spacer);
								MrBayesPrint ("%s                     selection is fixed to %1.2lf.\n", spacer, modelParams[i].ny98omega3Fixed);
								}
							}
						else if (!strcmp(modelParams[i].omegaVar, "M3"))
							{
							if (!strcmp(modelParams[i].m3omegapr, "Exponential"))
								{
								MrBayesPrint ("%s                     Nonsynonymous and synonymous rates for the tree classes of\n", spacer);
								MrBayesPrint ("%s                     omega are exponentially distributed random variables.\n", spacer);
								}
							else if (!strcmp(modelParams[i].m3omegapr, "Fixed"))
								{
								MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for the three omega\n", spacer);
								MrBayesPrint ("%s                     are fixed to %1.2lf, %1.2lf, and %1.2lf.\n", spacer, modelParams[i].m3omegaFixed[0], modelParams[i].m3omegaFixed[1], modelParams[i].m3omegaFixed[2]);
								}
							}
						else if (!strcmp(modelParams[i].omegaVar, "M10"))
							{
							MrBayesPrint ("%s                     Nonsynonymous/synonymous rate ratio for purifying \n", spacer);
							MrBayesPrint ("%s                     selection (class 1) has a Beta(alpha1,beta1) on the \n", spacer);
							MrBayesPrint ("%s                     interval (0,1). Nonsynonymous/synonymous rate ratio \n", spacer);
							MrBayesPrint ("%s                     for positive selection (class 2) has an offset \n", spacer);
							MrBayesPrint ("%s                     Gamma(alpha2,beta2) on the interval (1,Infinity).\n", spacer);
							}
							
						/* genetic code that is used (if nucmodel=codon) */
						MrBayesPrint ("%s         Code      = %s\n", spacer, modelParams[i].geneticCode);
						
						}
						
					}
				/* amino acid characters in this partition */
				else if (modelParams[i].dataType == PROTEIN)
					{
					/* constraints on rates of substitution in 20 X 20 matrix */
					if (!strcmp(modelParams[i].aaModelPr, "Mixed"))
						MrBayesPrint ("%s         Aamodel   = Mixture of models with fixed rate matrices\n", spacer);
					else
						MrBayesPrint ("%s         Aamodel   = %s\n", spacer, modelParams[i].aaModel);
					/* revmat rates */
					if (!strcmp(modelParams[i].aaModelPr, "Fixed") && !strcmp(modelParams[i].aaModel, "Gtr"))
						{
						if (!strcmp(modelParams[i].aaRevMatPr,"Dirichlet"))
							{
							for (j=0; j<190; j++)
								if (AreDoublesEqual(modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[j], 0.00001) == NO)
									break;
							if (j == 190)
								{
								MrBayesPrint ("%s                     Substitution rates have a Dirichlet(%1.2lf,%1.2lf,...) prior\n",
									spacer, modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[0]);
								}
							else
								{
								MrBayesPrint ("%s                     Substitution rates have a Dirichlet(\n", spacer);
								for (j=0; j<190; j++)
									{
									if (j % 10 == 0)
										MrBayesPrint ("%s                        ", spacer);
									MrBayesPrint ("%1.2lf", modelParams[i].aaRevMatDir[j]);
									if (j == 189)
										MrBayesPrint (") prior\n");
									else if ((j+1) % 10 == 0)
										MrBayesPrint (",\n");
									else
										MrBayesPrint (",");
									}
								}
							}
						else /* if (!strcmp(modelParams[i].aaRevMatPr,"Fixed")) */
							{
							for (j=0; j<190; j++)
								if (AreDoublesEqual(modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[j], 0.00001) == NO)
									break;
							if (j == 190)
								{
								MrBayesPrint ("%s                     Substitution rates are fixed to (%1.1lf,%1.1lf,...)\n",
									spacer, modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[0]);
								}
							else
								{
								MrBayesPrint ("%s                     Substitution rates are fixed to (\n", spacer);
								for (j=0; j<190; j++)
									{
									if (j % 10 == 0)
										MrBayesPrint ("%s                        ", spacer);
									MrBayesPrint ("%1.1lf", modelParams[i].aaRevMatFix[j]);
									if (j == 189)
										MrBayesPrint (") prior\n");
									else if ((j+1) % 10 == 0)
										MrBayesPrint (",\n");
									else
										MrBayesPrint (",");
									}
								}
							}
						}
					}
				/* restriction site or morphological characters in this partition */
				else if (modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD)
					{
					/* what type of characters are sampled? */
					MrBayesPrint ("%s         Coding    = %s\n", spacer, modelParams[i].coding);
					}
					
				/* is there rate variation in a single site across the tree? */
				if (((modelParams[i].dataType == DNA || modelParams[i].dataType == RNA) && !strcmp(modelParams[i].nucModel, "4by4")) || modelParams[i].dataType == PROTEIN)
					{
					/* do rates change on tree accoding to covarion model? */
					MrBayesPrint ("%s         Covarion  = %s\n", spacer, modelParams[i].covarionModel);
					if (!strcmp(modelParams[i].covarionModel, "Yes"))
						{
						/* distribution on switching parameters, if appropriate */
						if (!strcmp(modelParams[i].covSwitchPr,"Uniform"))
							{
							MrBayesPrint ("%s                     Switching rates have independent uniform dist-\n", spacer);
							MrBayesPrint ("%s                     ributions on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1]);
							}
						else if (!strcmp(modelParams[i].covSwitchPr,"Exponential"))
							{
							MrBayesPrint ("%s                     Switching rates have independent exponential\n", spacer);
							MrBayesPrint ("%s                     distributions with parameters (%1.2lf).\n", spacer, modelParams[i].covswitchExp);
							}
						else
							{
							MrBayesPrint ("%s                     Switching rates are fixed to %1.2lf and %1.2lf.\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[0]);
							}
						ns *= 2;
						}
					}

				/* now, let's deal with variation in omega */
				if ((modelParams[i].dataType == DNA || modelParams[i].dataType == RNA) && !strcmp(modelParams[i].nucModel,"Codon"))
					{
					MrBayesPrint ("%s         Omegavar  = %s\n", spacer, modelParams[i].omegaVar);
					if (!strcmp(modelParams[i].geneticCode, "Universal"))
						ns = 61;
					else if (!strcmp(modelParams[i].geneticCode, "Vertmt"))
						ns = 60;
					else if (!strcmp(modelParams[i].geneticCode, "Mycoplasma"))
						ns = 62;
					else if (!strcmp(modelParams[i].geneticCode, "Yeast"))
						ns = 62;
					else if (!strcmp(modelParams[i].geneticCode, "Ciliates"))
						ns = 63;
					else if (!strcmp(modelParams[i].geneticCode, "Metmt"))
						ns = 62;
					}

				/* what assumptions are made about the state frequencies? */
				if (modelParams[i].dataType != CONTINUOUS)
					{
					if (modelParams[i].dataType == STANDARD)
						MrBayesPrint ("%s         # States  = Variable, up to 10\n", spacer);
					else if (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA)
						{
						if (!strcmp(modelParams[i].nucModel, "Doublet"))
							MrBayesPrint ("%s         # States  = %d\n", spacer, 16);
						else if (!strcmp(modelParams[i].nucModel, "Codon"))	
							MrBayesPrint ("%s         # States  = %d\n", spacer, 61);
						else
							MrBayesPrint ("%s         # States  = %d\n", spacer, 4);
						}
					else
						MrBayesPrint ("%s         # States  = %d\n", spacer, ns);
					if (modelParams[i].dataType == STANDARD)
						{
						if (!strcmp(modelParams[i].symPiPr,"Fixed"))
						        { 
							  if (AreDoublesEqual(modelParams[i].symBetaFix, -1.0,ETA)==YES)
								MrBayesPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
							else
								MrBayesPrint ("%s                     Symmetric Dirichlet is fixed to %1.2lf\n", spacer, modelParams[i].symBetaFix);
							}
						else if (!strcmp(modelParams[i].symPiPr,"Uniform"))
							{
							MrBayesPrint ("%s                     Symmetric Dirichlet has a Uniform(%1.2lf,%1.2lf) prior\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1]);
							}
						else
							{
							MrBayesPrint ("%s                     Symmetric Dirichlet has a Exponential(%1.2lf) prior\n", spacer, modelParams[i].symBetaExp);
							}
						}
					else if (modelParams[i].dataType == RESTRICTION)
						{
						/* distribution on state frequencies for restriction site model */
						if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
							{
							MrBayesPrint ("%s                     State frequencies have a Dirichlet (%1.2lf,%1.2lf) prior\n", spacer,
								modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1]);
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
							{
							MrBayesPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
							{
							MrBayesPrint ("%s                     State frequencies have been fixed by the user\n", spacer);
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
							{
							MrBayesPrint ("%s                     State frequencies have been fixed to the empirical frequencies in the data\n", spacer);
							}
						}
					else if (modelParams[i].dataType == PROTEIN)
						{
						/* distribution on state frequencies for aminoacid model */
						if (!strcmp(modelParams[i].aaModelPr, "Fixed") && (strcmp(modelParams[i].aaModel, "Equalin")==0 ||
							strcmp(modelParams[i].aaModel, "Gtr")==0))
							{
							if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
								{
								MrBayesPrint ("%s                     State frequencies have a Dirichlet prior\n", spacer);
								MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,", spacer,
									modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1], modelParams[i].stateFreqsDir[2],
									modelParams[i].stateFreqsDir[3], modelParams[i].stateFreqsDir[4]);
								MrBayesPrint ("%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,\n",
									modelParams[i].stateFreqsDir[5], modelParams[i].stateFreqsDir[6], modelParams[i].stateFreqsDir[7],
									modelParams[i].stateFreqsDir[8], modelParams[i].stateFreqsDir[9]);
								MrBayesPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,", spacer,
									modelParams[i].stateFreqsDir[10], modelParams[i].stateFreqsDir[11], modelParams[i].stateFreqsDir[12],
									modelParams[i].stateFreqsDir[13], modelParams[i].stateFreqsDir[14]);
								MrBayesPrint ("%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n",
									modelParams[i].stateFreqsDir[15], modelParams[i].stateFreqsDir[16], modelParams[i].stateFreqsDir[17],
									modelParams[i].stateFreqsDir[18], modelParams[i].stateFreqsDir[19]);
								}
							else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
								{
								MrBayesPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
								}
							else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
								{
								MrBayesPrint ("%s                     State frequencies have been fixed by the user\n", spacer);
								}
							else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
								{
								MrBayesPrint ("%s                     State frequencies have been fixed to the empirical frequencies in the data\n", spacer);
								}
							}
						else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && !strcmp(modelParams[i].aaModel, "Poisson"))
							{
							MrBayesPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
							}
						else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && strcmp(modelParams[i].aaModel, "Equalin") && strcmp(modelParams[i].aaModel, "Poisson"))
							{
							MrBayesPrint ("%s                     State frequencies are fixed to the %s frequencies\n", spacer, modelParams[i].aaModel);
							}
						else
							{
							MrBayesPrint ("%s                     State frequencies come from the mixture of models\n", spacer);
							}
						}
					else
						{
						/* distribution on state frequencies for all other models */
						if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
							{
							MrBayesPrint ("%s                     State frequencies have a Dirichlet prior\n", spacer);
							if (!strcmp(modelParams[i].nucModel, "Doublet"))
								{
								MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf,\n", spacer,
									modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1], modelParams[i].stateFreqsDir[2],
									modelParams[i].stateFreqsDir[3]);
								MrBayesPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf,\n", spacer,
									modelParams[i].stateFreqsDir[4], modelParams[i].stateFreqsDir[5], modelParams[i].stateFreqsDir[6],
									modelParams[i].stateFreqsDir[7]);
								MrBayesPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf,\n", spacer,
									modelParams[i].stateFreqsDir[8], modelParams[i].stateFreqsDir[9], modelParams[i].stateFreqsDir[10],
									modelParams[i].stateFreqsDir[11]);
								MrBayesPrint ("%s                     %1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
									modelParams[i].stateFreqsDir[12], modelParams[i].stateFreqsDir[13], modelParams[i].stateFreqsDir[14],
									modelParams[i].stateFreqsDir[15]);
								}
							else if (!strcmp(modelParams[i].nucModel, "4by4"))
								{
								MrBayesPrint ("%s                     (%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer,
									modelParams[i].stateFreqsDir[0], modelParams[i].stateFreqsDir[1], modelParams[i].stateFreqsDir[2],
									modelParams[i].stateFreqsDir[3]);
								}
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
							{
							MrBayesPrint ("%s                     State frequencies are fixed to be equal\n", spacer);
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
							{
							MrBayesPrint ("%s                     State frequencies have been fixed by the user\n", spacer);
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
							{
							MrBayesPrint ("%s                     State frequencies have been fixed to the empirical frequencies in the data\n", spacer);
							}
						}
					}
				else
					MrBayesPrint ("%s         # States  = Infinity\n", spacer);

				/* now, let's deal with rate variation across sites */
				if (modelParams[i].dataType != CONTINUOUS)
					{
					if (((modelParams[i].dataType == DNA || modelParams[i].dataType == RNA) && strcmp(modelParams[i].nucModel,"Codon")) ||
					     modelParams[i].dataType == PROTEIN || modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD)
						{
						if (!strcmp(modelParams[i].covarionModel, "No"))
							MrBayesPrint ("%s         Rates     = %s\n", spacer, modelParams[i].ratesModel);
						else
							{
							if (!strcmp(modelParams[i].ratesModel, "Propinv"))
								MrBayesPrint ("%s         Rates     = Equal ", spacer);
							else if (!strcmp(modelParams[i].ratesModel, "Invgamma"))
								MrBayesPrint ("%s         Rates     = Gamma ", spacer);
							else
								MrBayesPrint ("%s         Rates     = %s ", spacer, modelParams[i].ratesModel);
							MrBayesPrint ("(+ Propinv induced by covarion model)\n");
							}
						
						if ((modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD) && !strcmp(modelParams[i].ratesModel, "Adgamma"))
							{
							
							}
						else
							{
							if ((!strcmp(modelParams[i].ratesModel, "Invgamma") || !strcmp(modelParams[i].ratesModel, "Gamma") || !strcmp(modelParams[i].ratesModel, "Adgamma")))
								{
								/* distribution on shape parameter, if appropriate */
								if (!strcmp(modelParams[i].shapePr,"Uniform"))
									{
									MrBayesPrint ("%s                     Gamma shape parameter is uniformly dist-\n", spacer);
									MrBayesPrint ("%s                     ributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1]);
									}
								else if (!strcmp(modelParams[i].shapePr,"Exponential"))
									{
									MrBayesPrint ("%s                     Gamma shape parameter is exponentially\n", spacer);
									MrBayesPrint ("%s                     distributed with parameter (%1.2lf).\n", spacer, modelParams[i].shapeExp);
									}
								else
									{
									MrBayesPrint ("%s                     Gamma shape parameter is fixed to %1.2lf.\n", spacer, modelParams[i].shapeFix);
									}
								}
							if ((!strcmp(modelParams[i].ratesModel, "Propinv") || !strcmp(modelParams[i].ratesModel, "Invgamma")) && !strcmp(modelParams[i].covarionModel, "No"))
								{
								/* distribution on pInvar parameter, if appropriate */
								if (!strcmp(modelParams[i].pInvarPr,"Uniform"))
									{
									MrBayesPrint ("%s                     Proportion of invariable sites is uniformly dist-\n", spacer);
									MrBayesPrint ("%s                     ributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1]);
									}
								else
									{
									MrBayesPrint ("%s                     Proportion of invariable sites is fixed to %1.2lf.\n", spacer, modelParams[i].pInvarFix);
									}
								}
							if (!strcmp(modelParams[i].ratesModel, "Adgamma"))
								{
								/* distribution on correlation parameter, if appropriate */
								if (!strcmp(modelParams[i].adGammaCorPr,"Uniform"))
									{
									MrBayesPrint ("%s                     Rate correlation parameter is uniformly dist-\n", spacer);
									MrBayesPrint ("%s                     ributed on the interval (%1.2lf,%1.2lf).\n", spacer, modelParams[i].corrUni[0], modelParams[i].corrUni[1]);
									}
								else
									{
									MrBayesPrint ("%s                     Rate correlation parameter is fixed to %1.2lf.\n", spacer, modelParams[i].corrFix);
									}
								}
							
							if (!strcmp(modelParams[i].ratesModel, "Gamma") || !strcmp(modelParams[i].ratesModel, "Invgamma") || !strcmp(modelParams[i].ratesModel, "Adgamma"))
								{
								/* how many categories is the continuous gamma approximated by? */
								MrBayesPrint ("%s                     Gamma distribution is approximated using %d categories.\n", spacer, modelParams[i].numGammaCats);
								}
							}
						
						}
					}

				}
			/* end description of discrete models */
			}

		if (i != numCurrentDivisions - 1)
			MrBayesPrint ("\n");
		
		}

#	if 1
	MrBayesPrint ("\n");
	MrBayesPrint ("%s   Active parameters: \n\n", spacer);
	if (numCurrentDivisions > 1)
		{ 
		MrBayesPrint ("%s                       Partition(s)\n", spacer);
		MrBayesPrint ("%s      Parameters     ", spacer);
		for (i=0; i<numCurrentDivisions; i++)
			MrBayesPrint (" %2d", i+1);
		MrBayesPrint ("\n");
		MrBayesPrint ("%s      ---------------", spacer);
		for (i=0; i<numCurrentDivisions; i++)
			MrBayesPrint ("---");
		MrBayesPrint ("\n");
		}
	else
		{
		MrBayesPrint ("%s      Parameters\n", spacer);
		MrBayesPrint ("%s      ------------------\n", spacer);
		}
	for (j=0; j<NUM_LINKED; j++)
		{
		
		shouldPrint = NO;
		for (i=0; i<numCurrentDivisions; i++)
			{
			if (activeParams[j][i] == -1)
				{}
			else
				shouldPrint = YES;
			}
		if (shouldPrint == NO)
			continue;
		
		if (j == P_TRATIO)
			{
			MrBayesPrint ("%s      Tratio         ", spacer);
			}
		else if (j == P_REVMAT)
			{
			MrBayesPrint ("%s      Revmat         ", spacer);
			}
		else if (j == P_OMEGA)
			{
			MrBayesPrint ("%s      Omega          ", spacer);
			}
		else if (j == P_PI)
			{
			MrBayesPrint ("%s      Statefreq      ", spacer);
			}
		else if (j == P_SHAPE)
			{
			MrBayesPrint ("%s      Shape          ", spacer);
			}
		else if (j == P_PINVAR)
			{
			MrBayesPrint ("%s      Pinvar         ", spacer);
			}
		else if (j == P_CORREL)
			{
			MrBayesPrint ("%s      Correlation    ", spacer);
			}
		else if (j == P_SWITCH)
			{
			MrBayesPrint ("%s      Switchrates    ", spacer);
			}
		else if (j == P_RATEMULT)
			{
			MrBayesPrint ("%s      Ratemultiplier ", spacer);
			}
		else if (j == P_TOPOLOGY)
			{
			MrBayesPrint ("%s      Topology       ", spacer);
			}
		else if (j == P_BRLENS)
			{
			MrBayesPrint ("%s      Brlens         ", spacer);
			}
		else if (j == P_SPECRATE)
			{
			MrBayesPrint ("%s      Speciationrate ", spacer);
			}
		else if (j == P_EXTRATE)
			{
			MrBayesPrint ("%s      Extinctionrate ", spacer);
			}
		else if (j == P_THETA)
			{
			MrBayesPrint ("%s      Theta          ", spacer);
			}
		else if (j == P_GENEMU) /*prior for joint gene tree */
			{
			MrBayesPrint ("%s      GeneMu          ", spacer);
			}
		else if (j == P_GROWTH)
			{
			MrBayesPrint ("%s      Growthrate     ", spacer);
			} 
		else if (j == P_AAMODEL)
			{
			MrBayesPrint ("%s      Aamodel        ", spacer);
			}
		else if (j == P_BRCORR)
			{
			MrBayesPrint ("%s      Brownian corr. ", spacer);
			}
		else if (j == P_BRSIGMA)
			{
			MrBayesPrint ("%s      Brownian sigma ", spacer);
			}
			
		for (i=0; i<numCurrentDivisions; i++)
			{
			if (activeParams[j][i] == -1)
				MrBayesPrint ("  .");
			else
				MrBayesPrint (" %2d", activeParams[j][i]);
			}
		MrBayesPrint ("\n");
		}
	if (numCurrentDivisions > 1)
		{ 
		MrBayesPrint ("%s      ---------------", spacer);
		for (i=0; i<numCurrentDivisions; i++)
			MrBayesPrint ("---");
		MrBayesPrint ("\n");
		}
	else
		{
		MrBayesPrint ("%s      ------------------\n", spacer);
		}

	MrBayesPrint ("\n");

	if (numCurrentDivisions > 1)
		MrBayesPrint ("%s      All parameters can be linked or unlinked across partitions\n\n", spacer);
		
	nOfParam = 0;
	for (j=0; j<NUM_LINKED; j++)
		for (i=0; i<numCurrentDivisions; i++)
			if (activeParams[j][i] > nOfParam)
				nOfParam = activeParams[j][i];
	nOfParam++;

	for (k=0; k<nOfParam; k++)
		{
		isFirst = YES;
		for (j=0; j<NUM_LINKED; j++)
			{
			for (i=0; i<numCurrentDivisions; i++)
				{
				if (activeParams[j][i] == k && isFirst == YES)
					{
					isFirst = NO;

					MrBayesPrint ("%s   %4d -- ", spacer, k);
					
					if (j == P_TRATIO)
						MrBayesPrint (" Parameter  = Tratio\n", spacer);
					else if (j == P_REVMAT)
						MrBayesPrint (" Parameter  = Revmat\n", spacer);
					else if (j == P_OMEGA)
						MrBayesPrint (" Parameter  = Omega\n", spacer);
					else if (j == P_PI)
						MrBayesPrint (" Parameter  = Statefreq\n", spacer);
					else if (j == P_SHAPE)
						MrBayesPrint (" Parameter  = Shape\n", spacer);
					else if (j == P_PINVAR)
						MrBayesPrint (" Parameter  = Pinvar\n", spacer);
					else if (j == P_CORREL)
						MrBayesPrint (" Parameter  = Correlation\n", spacer);
					else if (j == P_SWITCH)
						MrBayesPrint (" Parameter  = Switchrates\n", spacer);
					else if (j == P_RATEMULT)
						MrBayesPrint (" Parameter  = Ratemultiplier\n", spacer);
					else if (j == P_TOPOLOGY)
						MrBayesPrint (" Parameter  = Topology\n", spacer);
					else if (j == P_BRLENS)
						MrBayesPrint (" Parameter  = Brlens\n", spacer);
					else if (j == P_SPECRATE)
						MrBayesPrint (" Parameter  = Speciationrate\n", spacer);
					else if (j == P_EXTRATE)
						MrBayesPrint (" Parameter  = Extinctionrate\n", spacer);
					else if (j == P_THETA)
						{
						if (!strcmp(modelParams[i].ploidy,"Haploid"))
							MrBayesPrint (" Parameter  = Theta (2Nu)\n", spacer);
						else
							MrBayesPrint (" Parameter  = Theta (4Nu)\n", spacer);
						}
					else if (j == P_GENEMU) /*prior for joint gene tree */
                        			MrBayesPrint (" Parameter  = GeneMu\n", spacer); 
					else if (j == P_GROWTH)
						MrBayesPrint (" Parameter  = Growthrate (r/u)\n", spacer); 
					else if (j == P_AAMODEL)
						MrBayesPrint (" Parameter  = Aamodel\n", spacer);
					else if (j == P_BRCORR)
						MrBayesPrint (" Parameter  = Browncorr\n", spacer);
					else if (j == P_BRSIGMA)
						MrBayesPrint (" Parameter  = Brownscale\n", spacer);

					if (j == P_TRATIO)
						{
						if (!strcmp(modelParams[i].tRatioPr,"Beta"))
							MrBayesPrint ("%s            Prior      = Beta(%1.2lf,%1.2lf)\n", spacer, modelParams[i].tRatioDir[0], modelParams[i].tRatioDir[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].tRatioFix);
						}
					else if (j == P_REVMAT)
						{
						if (!strcmp(modelParams[i].revMatPr,"Dirichlet"))
							MrBayesPrint ("%s            Prior      = Dirichlet(%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer, 
							modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
							modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer, 
							modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2],
							modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);
						}
					else if (j == P_OMEGA)
						{
						if (!strcmp(modelParams[i].omegaVar,"Equal"))
							{
							if (!strcmp(modelParams[i].omegaPr,"Dirichlet"))
								MrBayesPrint ("%s            Prior      = Dirichlet(%1.2lf,%1.2lf)\n", spacer, modelParams[i].omegaDir[0], modelParams[i].omegaDir[1]);
							else
								MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].omegaFix);
							}
						else if (!strcmp(modelParams[i].omegaVar,"Ny98"))
							{
							if (!strcmp(modelParams[i].ny98omega1pr,"Beta"))
								MrBayesPrint ("%s            Prior      = Beta(%1.2lf,%1.2lf) on omega(1)\n", spacer, modelParams[i].ny98omega1Beta[0], modelParams[i].ny98omega1Beta[1]);
							else
								MrBayesPrint ("%s            Prior      = Fixed(%1.2lf) on omega(1)\n", spacer, modelParams[i].ny98omega1Fixed);
							MrBayesPrint ("%s                         Fixed(1.00) on omega(2)\n", spacer);
							if (!strcmp(modelParams[i].ny98omega3pr,"Uniform"))
								MrBayesPrint ("%s                         Uniform(%1.2lf,%1.2lf) on omega(3)\n", spacer, modelParams[i].ny98omega3Uni[0], modelParams[i].ny98omega3Uni[1]);
							else if (!strcmp(modelParams[i].ny98omega3pr,"Exponential"))
								MrBayesPrint ("%s                         Exponential(%1.2lf) on omega(3)\n", spacer, modelParams[i].ny98omega3Exp);
							else
								MrBayesPrint ("%s                         Fixed(%1.2lf) on omega(3)\n", spacer, modelParams[i].ny98omega3Fixed);
							if (!strcmp(modelParams[i].codonCatFreqPr,"Dirichlet"))
								MrBayesPrint ("%s                         Dirichlet(%1.2lf,%1.2lf,%1.2lf) on pi(-), pi(N), and pi(+)\n", spacer, modelParams[i].codonCatDir[0], modelParams[i].codonCatDir[1], modelParams[i].codonCatDir[2]);
							else
								MrBayesPrint ("%s                         Fixed(%1.2lf,%1.2lf,%1.2lf) on pi(-), pi(N), and pi(+)\n", spacer, modelParams[i].codonCatFreqFix[0], modelParams[i].codonCatFreqFix[1], modelParams[i].codonCatFreqFix[2]);
							}
						else if (!strcmp(modelParams[i].omegaVar,"M3"))
							{
							if (!strcmp(modelParams[i].m3omegapr,"Exponential"))
								MrBayesPrint ("%s                         dN1, dN2, dN3, and dS are all exponential random variables\n", spacer);
							else
								MrBayesPrint ("%s                         Fixed(%1.2lf,%1.2lf,%1.2lf) for three selection categories\n", spacer, modelParams[i].m3omegaFixed[0], modelParams[i].m3omegaFixed[1], modelParams[i].m3omegaFixed[2]);
							if (!strcmp(modelParams[i].codonCatFreqPr,"Dirichlet"))
								MrBayesPrint ("%s                         Dirichlet(%1.2lf,%1.2lf,%1.2lf) on pi(1), pi(2), and pi(3)\n", spacer, modelParams[i].codonCatDir[0], modelParams[i].codonCatDir[1], modelParams[i].codonCatDir[2]);
							else
								MrBayesPrint ("%s                         Fixed(%1.2lf,%1.2lf,%1.2lf) on pi(1), pi(2), and pi(3)\n", spacer, modelParams[i].codonCatFreqFix[0], modelParams[i].codonCatFreqFix[1], modelParams[i].codonCatFreqFix[2]);
							}
						else if (!strcmp(modelParams[i].omegaVar,"M10"))
							{
							MrBayesPrint ("%s                         With probability pi(1), omega is drawn from a Beta(alpha1,beta1) \n", spacer);
							MrBayesPrint ("%s                         distribution and with probability pi(2), omega is drawn from\n", spacer);
							MrBayesPrint ("%s                         a Gamma(alpha2,beta2) distribution.\n", spacer);
							if (!strcmp(modelParams[i].m10betapr,"Uniform"))
								{
								MrBayesPrint ("%s                         The parameters of the beta distribution each follow \n", spacer);
								MrBayesPrint ("%s                         independent Uniform(%1.2lf,%1.2lf) priors\n", spacer, modelParams[i].m10betaUni[0], modelParams[i].m10betaUni[1]);
								}
							else if (!strcmp(modelParams[i].m10betapr,"Exponential"))
								{
								MrBayesPrint ("%s                         The parameters of the beta distribution each follow \n", spacer);
								MrBayesPrint ("%s                         independent Exponential(%1.2lf) priors\n", spacer, modelParams[i].m10betaExp);
								}
							else
								{
								MrBayesPrint ("%s                         The parameters of the beta distribution are fixed to \n", spacer);
								MrBayesPrint ("%s                         %1.2lf and %1.2lf\n", spacer, modelParams[i].m10betaFix[0], modelParams[i].m10betaFix[0]);
								}

							if (!strcmp(modelParams[i].m10gammapr,"Uniform"))
								{
								MrBayesPrint ("%s                         The parameters of the gamma distribution each follow  \n", spacer);
								MrBayesPrint ("%s                         independent Uniform(%1.2lf,%1.2lf) priors\n", spacer, modelParams[i].m10gammaUni[0], modelParams[i].m10gammaUni[1]);
								}
							else if (!strcmp(modelParams[i].m10gammapr,"Exponential"))
								{
								MrBayesPrint ("%s                         The parameters of the gamma distribution each follow \n", spacer);
								MrBayesPrint ("%s                         independent Exponential(%1.2lf) priors\n", spacer, modelParams[i].m10gammaExp);
								}
							else
								{
								MrBayesPrint ("%s                         The parameters of the gamma distribution are fixed to \n", spacer);
								MrBayesPrint ("%s                         %1.2lf and %1.2lf\n", spacer, modelParams[i].m10gammaFix[0], modelParams[i].m10gammaFix[0]);
								}

							if (!strcmp(modelParams[i].codonCatFreqPr,"Dirichlet"))
								MrBayesPrint ("%s                         Dirichlet(%1.2lf,%1.2lf) on pi(1) and pi(2)\n", spacer, modelParams[i].codonCatDir[0], modelParams[i].codonCatDir[1]);
							else
								MrBayesPrint ("%s                         Fixed(%1.2lf,%1.2lf) on pi(1) and pi(2)\n", spacer, modelParams[i].codonCatFreqFix[0], modelParams[i].codonCatFreqFix[1]);
							}
						}
					else if (j == P_PI)
						{
						if (modelParams[i].dataType == STANDARD)
							{
							if (!strcmp(modelParams[i].symPiPr, "Uniform"))
								MrBayesPrint ("%s            Prior      = Symmetric dirichlet with uniform(%1.2lf,%1.2lf) variance parameter\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1]);
							else if (!strcmp(modelParams[i].symPiPr, "Exponential"))
								MrBayesPrint ("%s            Prior      = Symmetric dirichlet with exponential(%1.2lf) variance parameter\n", spacer, modelParams[i].symBetaExp);
							else
							        { /* modelParams[i].symBetaFix == -1 */
								  if (AreDoublesEqual(modelParams[i].symBetaFix, 1.0, ETA)==YES)
									MrBayesPrint ("%s            Prior      = State frequencies are equal\n", spacer);
								else
									MrBayesPrint ("%s            Prior      = Symmetric dirichlet with fixed(%1.2lf) variance parameter\n", spacer, modelParams[i].symBetaFix);
								}
							}
						else if (modelParams[i].dataType == PROTEIN)
							{
							if (!strcmp(modelParams[i].aaModelPr, "Fixed") && !strcmp(modelParams[i].aaModel, "Equalin"))
								{
								if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
									{
									MrBayesPrint ("%s            Prior      = Dirichlet\n", spacer);
									}
								else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Equal"))
									{
									MrBayesPrint ("%s            Prior      = Fixed (equal frequencies)\n", spacer);
									}
								else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"User"))
									{
									MrBayesPrint ("%s            Prior      = Fixed (user-specified)\n", spacer);
									}
								else if (!strcmp(modelParams[i].stateFreqPr,"Fixed") && !strcmp(modelParams[i].stateFreqsFixType,"Empirical"))
									{
									MrBayesPrint ("%s            Prior      = Fixed (empirical frequencies)\n", spacer);
									}
								}
							else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && !strcmp(modelParams[i].aaModel, "Poisson"))
								{
								MrBayesPrint ("%s            Prior      = Fixed (equal frequencies)\n", spacer);
								}
							else if (!strcmp(modelParams[i].aaModelPr, "Fixed") && strcmp(modelParams[i].aaModel, "Equalin") && strcmp(modelParams[i].aaModel, "Poisson"))
								{
								MrBayesPrint ("%s            Prior      = Fixed (%s frequencies)\n", spacer, modelParams[i].aaModel);
								}
							else
								{
								MrBayesPrint ("%s            Prior      = Fixed (from mixture of models)\n", spacer);
								}
							}
						else
							{
							if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
								MrBayesPrint ("%s            Prior      = Dirichlet\n", spacer);
							else
								MrBayesPrint ("%s            Prior      = Fixed\n", spacer);
							}
						}
					else if (j == P_SHAPE)
						{
						if (!strcmp(modelParams[i].shapePr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1]);
						else if (!strcmp(modelParams[i].shapePr,"Exponential"))
							MrBayesPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, modelParams[i].shapeExp);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].shapeFix);
						}
					else if (j == P_PINVAR)
						{
						if (!strcmp(modelParams[i].pInvarPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].pInvarFix);
						}
					else if (j == P_CORREL)
						{
						if (!strcmp(modelParams[i].adGammaCorPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].corrUni[0], modelParams[i].corrUni[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].corrFix);
						}
					else if (j == P_SWITCH)
						{
						if (!strcmp(modelParams[i].covSwitchPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1]);
						else if (!strcmp(modelParams[i].covSwitchPr,"Exponential"))
							MrBayesPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, modelParams[i].covswitchExp);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf,%1.2lf)\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[1]);
						}
					else if (j == P_RATEMULT)
						{
						MrBayesPrint ("%s            Prior      = Dirichlet(", spacer);
						for (d=n=0; d<numCurrentDivisions; d++)
							{
							if (activeParams[j][d] == k)
								n++;
							}
						for (d=m=0; d<numCurrentDivisions; d++)
							{
							if (activeParams[j][d] == k)
								{
								m++;
								if (m < n)
									MrBayesPrint ("%1.2lf,", modelParams[d].ratePrDir);
								else
									MrBayesPrint ("%1.2lf)\n", modelParams[d].ratePrDir);
								}
							}
						}
					else if (j == P_TOPOLOGY)
						{
						if (!strcmp(modelParams[i].topologyPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = All topologies equally probable a priori\n", spacer);
						else
							MrBayesPrint ("%s            Prior      = Prior on topologies obeys constraints\n", spacer);
						}
					else if (j == P_BRLENS)
						{
						if (!strcmp(modelParams[i].parsModel, "Yes"))
							MrBayesPrint ("%s            Prior      = Branch lengths reconstructed using parsimony\n", spacer);
						else
							{
							if (!strcmp(modelParams[i].brlensPr,"Unconstrained"))
								{
								MrBayesPrint ("%s            Prior      = Branch lengths are Unconstrained:%s", spacer, modelParams[i].unconstrainedPr);
								if (!strcmp(modelParams[i].unconstrainedPr, "Uniform"))
									MrBayesPrint ("(%1.1lf,%1.1lf)\n", modelParams[i].brlensUni[0], modelParams[i].brlensUni[1]);
								else
									MrBayesPrint ("(%1.1lf)\n", modelParams[i].brlensExp);
								}
							else
								{
								MrBayesPrint ("%s            Prior      = Branch lengths are Clock:%s\n", spacer, modelParams[i].clockPr);
								if (strcmp(modelParams[i].clockPr, "Coalescence"))
									{
									if (!strcmp(modelParams[i].treeHeightPr, "Gamma"))
										MrBayesPrint ("%s                         Tree height is Gamma(%1.3lf,%1.3lf) distributed\n", spacer, modelParams[i].treeHeightGamma[0], modelParams[i].treeHeightGamma[1]);
									else
										MrBayesPrint ("%s                         Tree height has an Exponential(%1.3lf) distribution\n", spacer, modelParams[i].treeHeightExp);
									}
								b = 0;
								for (a=0; a<numTaxa; a++)
									{
									if (taxaInfo[a].isDeleted == NO && taxaInfo[a].calibration.prior != unconstrained)
										b++;
									}
								for (a=0; a<MAX_NUM_CONSTRAINTS; a++)
									{
									if (modelParams[i].activeConstraints[a] == YES && constraintCalibration[a].prior != unconstrained)
										b++;
									}
								if (b > 0)
									{
									MrBayesPrint ("%s                         Branch lengths are constrained by the following age constraint%s:\n", spacer, b > 1 ? "s" : "");
									for (a=0; a<numTaxa; a++)
										{
										if (taxaInfo[a].isDeleted == NO && taxaInfo[a].calibration.prior != unconstrained)
											{
											GetNameFromString (taxaNames, tempName, a+1);
											MrBayesPrint ("%s                         -- The age of terminal \"%s\" is %s\n", spacer, tempName,
												taxaInfo[a].calibration.name);
											}
										}
									for (a=0; a<numDefinedConstraints; a++)
										{
										if (modelParams[i].activeConstraints[a-numTaxa] == YES && constraintCalibration[a].prior != unconstrained)
											{
											GetNameFromString (constraintNames, tempName, a+1);
											MrBayesPrint ("%s                         -- The age of constrained node \"%s\" is %s\n", spacer, tempName,
												constraintCalibration[a].name);
											}
										}
									}
								}
							}
						}
					else if (j == P_SPECRATE)
						{
						if (!strcmp(modelParams[i].speciationPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].speciationUni[0], modelParams[i].speciationUni[1]);
						else if (!strcmp(modelParams[i].speciationPr,"Exponential"))
							MrBayesPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, modelParams[i].speciationExp);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].speciationFix);
						}
					else if (j == P_EXTRATE)
						{
						if (!strcmp(modelParams[i].extinctionPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].extinctionUni[0], modelParams[i].extinctionUni[1]);
						else if (!strcmp(modelParams[i].extinctionPr,"Exponential"))
							MrBayesPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, modelParams[i].extinctionExp);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].extinctionFix);
						}
					else if (j == P_THETA)
						{
						if (!strcmp(modelParams[i].thetaPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].thetaUni[0], modelParams[i].thetaUni[1]);
						else if (!strcmp(modelParams[i].thetaPr,"Exponential"))
							MrBayesPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, modelParams[i].thetaExp);
						/*prior for joint gene tree*/
						else if (!strcmp(modelParams[i].thetaPr,"Gamma"))
                                                        MrBayesPrint ("%s            Prior      = Gamma(%1.4lf,%2.1lf)\n", spacer, modelParams[i].thetaGamma[0], modelParams[i].thetaGamma[1]);
						else if (!strcmp(modelParams[i].thetaPr,"Invgamma"))
                                                        MrBayesPrint ("%s            Prior      = InvGamma(%2.1lf,%1.4lf)\n", spacer, modelParams[i].thetaGamma[0], modelParams[i].thetaGamma[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].thetaFix);
						}
					else if (j == P_GENEMU)  /*prior for joint gene tree */
						{
						if (!strcmp(modelParams[i].GeneMuPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].GeneMuUni[0], modelParams[i].GeneMuUni[1]);
						else if (!strcmp(modelParams[i].GeneMuPr,"Gamma"))
							MrBayesPrint ("%s            Prior      = Gamma(%1.2lf,%1.2lf)\n", spacer, modelParams[i].GeneMuGamma[0],modelParams[i].GeneMuGamma[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].GeneMuFix);
						}
					else if (j == P_GROWTH)
						{
						if (!strcmp(modelParams[i].growthPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].growthUni[0], modelParams[i].growthUni[1]);
						else if (!strcmp(modelParams[i].growthPr,"Exponential"))
							MrBayesPrint ("%s            Prior      = Exponential(%1.2lf)\n", spacer, modelParams[i].growthExp);
						else if (!strcmp(modelParams[i].growthPr,"Normal"))
							MrBayesPrint ("%s            Prior      = Normal(%1.2lf,%1.2lf)\n", spacer, modelParams[i].growthNorm[0], modelParams[i].growthNorm[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].growthFix);
						} 
					else if (j == P_AAMODEL)
						{
						if (!strcmp(modelParams[i].aaModelPr,"Mixed"))
							{
							/* check to see if you have a uniform prior */
							isSame = YES;
							for (a=0; a<9; a++)
								for (b=a+1; b<10; b++)
								        /* modelParams[i].aaModelPrProbs[a] != modelParams[i].aaModelPrProbs[b] */
 						        if (AreDoublesEqual(modelParams[i].aaModelPrProbs[a], modelParams[i].aaModelPrProbs[b], ETA)==FALSE)
										isSame = NO;
							MrBayesPrint ("%s            Prior      = Poisson, Jones, Dayhoff, Mtrev, Mtmam, Wag, Rtrev,\n", spacer);
							if (isSame == YES)
								MrBayesPrint ("%s                         Cprev, Vt, and Blosum models have equal prior probability\n", spacer);
							else
								MrBayesPrint ("%s                         Cprev, Vt, and Blosum models have unequal prior probability\n", spacer);
							}
						else
							MrBayesPrint ("%s            Prior      = Fixed(%s)\n", spacer, modelParams[i].aaModel);
						}
					else if (j == P_BRCORR)
						{
						if (!strcmp(modelParams[i].brownCorPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownCorrUni[0], modelParams[i].brownCorrUni[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].brownCorrFix);
						}
					else if (j == P_BRSIGMA)
						{
						if (!strcmp(modelParams[i].brownScalesPr,"Uniform"))
							MrBayesPrint ("%s            Prior      = Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownScalesUni[0], modelParams[i].brownScalesUni[1]);
						else if (!strcmp(modelParams[i].brownScalesPr,"Gammamean"))
							MrBayesPrint ("%s            Prior      = Gamma Mean=<char. ave.> Var=%1.2lf\n", spacer, modelParams[i].brownScalesGammaMean);
						else if (!strcmp(modelParams[i].brownScalesPr,"Gamma"))
							MrBayesPrint ("%s            Prior      = Gamma Mean=%lf Var=%1.2lf\n", spacer, modelParams[i].brownScalesGamma[0], modelParams[i].brownScalesGamma[1]);
						else
							MrBayesPrint ("%s            Prior      = Fixed(%1.2lf)\n", spacer, modelParams[i].brownScalesFix);
						}
						
					if (numCurrentDivisions > 1)
						{
						n = 0;
						for (m=0; m<numCurrentDivisions; m++)
							if (activeParams[j][m] == k)
								n++;
						if (n == 1)
							MrBayesPrint ("%s            Partition  = ", spacer);
						else
							MrBayesPrint ("%s            Partitions = ", spacer);
						paramCount = 0;
						for (m=0; m<numCurrentDivisions; m++)
							{
							if (activeParams[j][m] == k)
								{
								paramCount++;
								if (n == 2)
									{
									if (paramCount == n)
										MrBayesPrint (" and ");
									}
								else if (n > 2)
									{
									if (paramCount == n)
										MrBayesPrint (", and ");
									else if (paramCount != 1)
										MrBayesPrint (", ");
									}
								MrBayesPrint ("%d", m+1);
								}
							}
						MrBayesPrint ("\n");
						}


					}
				}

			}
		}
		
#	endif

	return (NO_ERROR);
	
}





int DoLink (void)

{

	int			i, j, newLine;
	
	MrBayesPrint ("%s   Linking\n", spacer);
	
	/* update status of linkTable */
	for (j=0; j<NUM_LINKED; j++)
		{
		newLine = YES;
		for (i=0; i<numCurrentDivisions; i++)
			{
			if (tempLinkUnlink[j][i] == YES)
				{
				if (newLine == YES)
					{
					linkNum++;
					newLine = NO;
					}
				linkTable[j][i] = linkNum;
				}
			}
		}
	
#	if 0
	for (j=0; j<NUM_LINKED; j++)
		{
		MrBayesPrint ("%s   ", spacer);
		for (i=0; i<numCurrentDivisions; i++)
			MrBayesPrint ("%d", linkTable[j][i]);
		MrBayesPrint ("\n");
		}
#	endif

	/* reinitialize the temporary table */
	for (j=0; j<NUM_LINKED; j++)
		for (i=0; i<MAX_NUM_DIVS; i++)
			tempLinkUnlink[j][i] = NO;

	return (NO_ERROR);
	
}





int DoLinkParm (char *parmName, char *tkn)

{

	int			i, j, tempInt;

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
		return (ERROR);
		}
		
	if (inValidCommand == YES)
		{
		for (j=0; j<NUM_LINKED; j++)
			for (i=0; i<MAX_NUM_DIVS; i++)
				tempLinkUnlink[j][i] = NO;
		inValidCommand = NO;
		}

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else if (expecting == Expecting(EQUALSIGN))
		{
		expecting = Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(LEFTPAR))
		{
		/* initialize tempLinkUnlinkVec to no */
		for (i=0; i<MAX_NUM_DIVS; i++)
			tempLinkUnlinkVec[i] = NO;
		fromI = toJ = -1;
		foundDash = NO;
		expecting = Expecting(NUMBER) | Expecting(ALPHA);
		}
	else if (expecting == Expecting(RIGHTPAR))
		{
		if (fromI != -1)
			tempLinkUnlinkVec[fromI-1] = YES;
		/* now copy tempLinkUnlinkVec to appropriate row of tempLinkUnlink */
		if (!strcmp(parmName, "Tratio"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_TRATIO][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Revmat"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_REVMAT][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Omega"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_OMEGA][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Statefreq"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_PI][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Shape"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_SHAPE][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Pinvar"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_PINVAR][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Correlation"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_CORREL][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Ratemultiplier"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_RATEMULT][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Switchrates"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_SWITCH][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Topology"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_TOPOLOGY][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Brlens"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_BRLENS][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Speciationrate"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_SPECRATE][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Extinctionrate"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_EXTRATE][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Theta"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_THETA][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "GeneMuPr")) /* prior for joint gene tree */
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_GENEMU][i] = tempLinkUnlinkVec[i];
			}
		else if (!strcmp(parmName, "Growthrate"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_GROWTH][i] = tempLinkUnlinkVec[i];
			} 
		else if (!strcmp(parmName, "Aamodel"))
			{
			for (i=0; i<numCurrentDivisions; i++)
				tempLinkUnlink[P_AAMODEL][i] = tempLinkUnlinkVec[i];
			}

		else
			{
			MrBayesPrint ("%s   Couldn't find parameter %s to link\n", spacer, parmName);
			}
		
		expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(COMMA))
		{
		foundComma = YES;
		expecting = Expecting(NUMBER);
		}
	else if (expecting == Expecting(ALPHA))
		{
		if (IsSame ("All", tkn) == DIFFERENT)
			{
			MrBayesPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
			return (ERROR);
			}
		for (i=0; i<numCurrentDivisions; i++)
			tempLinkUnlinkVec[i] = YES;
		expecting  = Expecting(RIGHTPAR);
		}
	else if (expecting == Expecting(NUMBER))
		{
		sscanf (tkn, "%d", &tempInt);
		if (tempInt > numCurrentDivisions)
			{
			MrBayesPrint ("%s   Partition delimiter is too large\n", spacer);
			return (ERROR);
			}
		if (fromI == -1)
			fromI = tempInt;
		else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
			{
			toJ = tempInt;
			for (i=fromI-1; i<toJ; i++)
				tempLinkUnlinkVec[i] = YES;
			fromI = toJ = -1;
			foundDash = NO;
			}
		else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
			{
			tempLinkUnlinkVec[fromI-1] = YES;
			fromI = tempInt;
			foundComma = NO;
			}
		expecting  = Expecting(COMMA);
		expecting |= Expecting(DASH);
		expecting |= Expecting(RIGHTPAR);
		}
	else if (expecting == Expecting(DASH))
		{
		foundDash = YES;
		expecting = Expecting(NUMBER);
		}
	else
		return (ERROR);

	return (NO_ERROR);
	
}





int DoLset (void)

{

	int			i, nApplied, lastActive=0;
	
	nApplied = NumActiveParts ();
	for (i=numCurrentDivisions; i>=0; i--)
		{
		if (activeParts[i] == YES)
			{
			lastActive = i;
			break;
			}
		}
			
	/* MrBayesPrint ("\n"); */
	if (numCurrentDivisions == 1)
		MrBayesPrint ("%s   Successfully set likelihood model parameters\n", spacer);
	else 
		{
		if (nApplied == numCurrentDivisions || nApplied == 0)
			{
			MrBayesPrint ("%s   Successfully set likelihood model parameters to all\n", spacer);
			MrBayesPrint ("%s   applicable data partitions \n", spacer);
			}
		else
			{
			MrBayesPrint ("%s   Successfully set likelihood model parameters to\n", spacer);
			if (nApplied == 1)
				MrBayesPrint ("%s   partition", spacer);
			else
				MrBayesPrint ("%s   partitions", spacer);
			for (i=0; i<numCurrentDivisions; i++)
				{
				if (activeParts[i] == YES)
					{
					if (i == lastActive && nApplied > 1)
						MrBayesPrint (" and %d", i+1);
					else
						MrBayesPrint (" %d", i+1);
					if (nApplied > 2 && i != lastActive)
						MrBayesPrint (",");
					}
				}
			MrBayesPrint (" (if applicable)\n");
			}
		}

	return (NO_ERROR);
	
}





int DoLsetParm (char *parmName, char *tkn)

{

	int			i, j, tempInt, nApplied;
	char		tempStr[100];

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
		return (ERROR);
		}
	if (inValidCommand == YES)
		{
		for (i=0; i<MAX_NUM_DIVS; i++)
			activeParts[i] = NO;
		inValidCommand = NO;
		}

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		/* set Applyto (Applyto) *************************************************************/
		if (!strcmp(parmName, "Applyto"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(LEFTPAR);
			else if (expecting == Expecting(LEFTPAR))
				{
				for (i=0; i<MAX_NUM_DIVS; i++)
					activeParts[i] = NO;
				fromI = toJ = -1;
				foundDash = NO;
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				if (fromI != -1)
					activeParts[fromI-1] = YES;
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else if (expecting == Expecting(COMMA))
				{
				foundComma = YES;
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (IsSame ("All", tkn) == DIFFERENT)
					{
					MrBayesPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
					return (ERROR);
					}
				for (i=0; i<numCurrentDivisions; i++)
					activeParts[i] = YES;
				expecting  = Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt > numCurrentDivisions)
					{
					MrBayesPrint ("%s   Partition delimiter is too large\n", spacer);
					return (ERROR);
					}
				if (fromI == -1)
					fromI = tempInt;
				else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
					{
					toJ = tempInt;
					for (i=fromI-1; i<toJ; i++)
						activeParts[i] = YES;
					fromI = toJ = -1;
					foundDash = NO;
					}
				else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
					{
					activeParts[fromI-1] = YES;
					fromI = tempInt;
					foundComma = NO;
					}
					
				expecting  = Expecting(COMMA);
				expecting |= Expecting(DASH);
				expecting |= Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}
		/* set Nucmodel (nucModel) ************************************************************/
		else if (!strcmp(parmName, "Nucmodel"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					tempInt = NO;
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].nucModel, tempStr);
							modelParams[i].nStates = NumStates (i);
							
							/* set state frequencies back to default */
							strcpy(modelParams[i].stateFreqPr, "Dirichlet");
							strcpy(modelParams[i].stateFreqsFixType, "Equal");
							for (j=0; j<200; j++)
								{
								modelParams[i].stateFreqsFix[j] = 0.0;   
								modelParams[i].stateFreqsDir[j] = 1.0;
								}    
							tempInt = YES;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Nucmodel to %s\n", spacer, modelParams[i].nucModel);
							else
								MrBayesPrint ("%s   Setting Nucmodel to %s for partition %d\n", spacer, modelParams[i].nucModel, i+1);
							}
						}
					if (tempInt == YES)
						MrBayesPrint ("%s   Set state frequency prior to default\n", spacer);
					}
				else
					{
					MrBayesPrint ("%s   Invalid DNA substitution model\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Nst (nst) **********************************************************************/
		else if (!strcmp(parmName, "Nst"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
			else if (expecting == Expecting(NUMBER) || expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].nst, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Nst to %s\n", spacer, modelParams[i].nst);
							else
								MrBayesPrint ("%s   Setting Nst to %s for partition %d\n", spacer, modelParams[i].nst, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Nst argument\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Ncat (numGammaCats) ************************************************************/
		else if (!strcmp(parmName, "Ngammacat"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt >= 2 && tempInt < MAX_GAMMA_CATS)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
							{
							modelParams[i].numGammaCats = tempInt;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Ngammacat to %d\n", spacer, modelParams[i].numGammaCats);
							else
								MrBayesPrint ("%s   Setting Ngammacat to %d for partition %d\n", spacer, modelParams[i].numGammaCats, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Ngammacat argument (should be between 2 and %d)\n", spacer, MAX_GAMMA_CATS);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set NumM10GammaCats (numM10GammaCats) ************************************************************/
		else if (!strcmp(parmName, "NumM10GammaCats"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt >= 2 && tempInt < MAX_GAMMA_CATS)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
							{
							modelParams[i].numM10GammaCats = tempInt;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting NumM10GammaCats to %d\n", spacer, modelParams[i].numM10GammaCats);
							else
								MrBayesPrint ("%s   Setting NumM10GammaCats to %d for partition %d\n", spacer, modelParams[i].numM10GammaCats, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid NumM10GammaCats argument (should be between 2 and %d)\n", spacer, MAX_GAMMA_CATS);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set NumM10BetaCats (numM10BetaCats) ************************************************************/
		else if (!strcmp(parmName, "NumM10BetaCats"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt >= 2 && tempInt < MAX_GAMMA_CATS)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType != CONTINUOUS))
							{
							modelParams[i].numM10BetaCats = tempInt;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting NumM10BetaCats to %d\n", spacer, modelParams[i].numM10BetaCats);
							else
								MrBayesPrint ("%s   Setting NumM10BetaCats to %d for partition %d\n", spacer, modelParams[i].numM10BetaCats, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid NumM10GammaCats argument (should be between 2 and %d)\n", spacer, MAX_GAMMA_CATS);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Nbetacat (numBetaCats) *****************************************************/
		else if (!strcmp(parmName, "Nbetacat"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt >= 2 && tempInt < MAX_GAMMA_CATS)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
							{
							modelParams[i].numBetaCats = tempInt;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Nbetacat to %d\n", spacer, modelParams[i].numBetaCats);
							else
								MrBayesPrint ("%s   Setting Nbetacat to %d for partition %d\n", spacer, modelParams[i].numBetaCats, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Nbetacat argument (should be between 2 and %d)\n", spacer, MAX_GAMMA_CATS);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Aamodel (aaModel) **************************************************************/
		else if (!strcmp(parmName, "Aamodel"))
			{
			MrBayesPrint ("%s   Aamodel argument for lset deprecated.\n", spacer);
			MrBayesPrint ("%s   Use 'prset aamodelpr=fixed(<aamodel>)' instead.\n", spacer);
			return (ERROR);
			}
		/* set Parsmodel (useParsModel) *******************************************************/
		else if (!strcmp(parmName, "Parsmodel"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							{
							if (!strcmp(tempStr, "Yes"))
								strcpy(modelParams[i].parsModel, "Yes");
							else
								strcpy(modelParams[i].parsModel, "No");

							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Parsmodel to %s\n", spacer, modelParams[i].parsModel);
							else
								MrBayesPrint ("%s   Setting Parsmodel to %s for partition %d\n", spacer, modelParams[i].parsModel, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for using (so-called) parsimony model\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}			
		/* set Augment (augmentData) **********************************************************/
		else if (!strcmp(parmName, "Augment"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
							{
							if (!strcmp(tempStr, "Yes"))
								strcpy(modelParams[i].augmentData, "Yes");
							else
								strcpy(modelParams[i].augmentData, "No");

							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Augmentdata to %s\n", spacer, modelParams[i].augmentData);
							else
								MrBayesPrint ("%s   Setting Augmentdata to %s for partition %d\n", spacer, modelParams[i].augmentData, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for data augmentation\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}			
		/* set Omegavar (wVarModel) ***********************************************************/
		else if (!strcmp(parmName, "Omegavar"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].omegaVar, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Omegavar to %s\n", spacer, modelParams[i].omegaVar);
							else
								MrBayesPrint ("%s   Setting Omegavar to %s for partition %d\n", spacer, modelParams[i].omegaVar, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid omega variation argument\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Code (codeModel) ***************************************************************/
		else if (!strcmp(parmName, "Code"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].geneticCode, tempStr);
							SetCode (i);
							modelParams[i].nStates = NumStates (i);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Code to %s\n", spacer, modelParams[i].geneticCode);
							else
								MrBayesPrint ("%s   Setting Code to %s for partition %d\n", spacer, modelParams[i].geneticCode, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid genetic code argument\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Ploidy (ploidy) ***************************************************************/
		else if (!strcmp(parmName, "Ploidy"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].ploidy, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting ploidy level to %s\n", spacer, modelParams[i].ploidy);
							else
								MrBayesPrint ("%s   Setting ploidy level to %s for partition %d\n", spacer, modelParams[i].ploidy, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid ploidy level argument\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Rates (ratesModel) *************************************************************/
		else if (!strcmp(parmName, "Rates"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
							{
							if (!strcmp(tempStr, "Adgamma") && (modelParams[i].dataType != DNA && modelParams[i].dataType != RNA && modelParams[i].dataType != PROTEIN))
								{
								/* we won't apply an adgamma model to anything but DNA, RNA, or PROTEIN data */
								}
							else if ((!strcmp(tempStr, "Propinv") ||  !strcmp(tempStr, "Invgamma")) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
								{
								/* we will not apply pinvar to standard or restriction site data */
								}
							else
								{
								strcpy(modelParams[i].ratesModel, tempStr);
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Rates to %s\n", spacer, modelParams[i].ratesModel);
								else
									MrBayesPrint ("%s   Setting Rates to %s for partition %d\n", spacer, modelParams[i].ratesModel, i+1);
								}
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Rates argument\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Covarion (covarionModel) *******************************************************/
		else if (!strcmp(parmName, "Covarion"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
							{
							strcpy(modelParams[i].covarionModel, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Covarion to %s\n", spacer, modelParams[i].covarionModel);
							else
								MrBayesPrint ("%s   Setting Covarion to %s for partition %d\n", spacer, modelParams[i].covarionModel, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Rates argument\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Coding (missingType) ***********************************************************/
		else if (!strcmp(parmName, "Coding"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD))
							{
							strcpy(modelParams[i].coding, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Coding to %s\n", spacer, modelParams[i].coding);
							else
								MrBayesPrint ("%s   Setting Coding to %s for partition %d\n", spacer, modelParams[i].coding, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for missing patterns\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
				
				
				
				
		else
			return (ERROR);
		}

	return (NO_ERROR);

}





int DoPrset (void)

{

	int			i, nApplied, lastActive=0;

	nApplied = NumActiveParts ();
	for (i=numCurrentDivisions; i>=0; i--)
		{
		if (activeParts[i] == YES)
			{
			lastActive = i;
			break;
			}
		}
			
	if (numCurrentDivisions == 1)
		MrBayesPrint ("%s   Successfully set prior model parameters\n", spacer);
	else 
		{
		if (nApplied == numCurrentDivisions || nApplied == 0)
			{
			MrBayesPrint ("%s   Successfully set prior model parameters to all\n", spacer);
			MrBayesPrint ("%s   applicable data partitions \n", spacer);
			}
		else
			{
			MrBayesPrint ("%s   Successfully set prior model parameters to\n", spacer);
			if (nApplied == 1)
				MrBayesPrint ("%s   partition", spacer);
			else
				MrBayesPrint ("%s   partitions", spacer);
			for (i=0; i<numCurrentDivisions; i++)
				{
				if (activeParts[i] == YES)
					{
					if (i == lastActive && nApplied > 1)
						MrBayesPrint (" and %d", i+1);
					else
						MrBayesPrint (" %d", i+1);
					if (nApplied > 2 && i != lastActive)
						MrBayesPrint (",");
					}
				}
			MrBayesPrint (" (if applicable)\n");
			}
		}

	return (NO_ERROR);
	
}





int DoPrsetParm (char *parmName, char *tkn)

{

	int			i, j, k, tempInt, nApplied, howMany, ns;
	MrBFlt		tempD, sum;
	char		tempStr[100];

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
		return (ERROR);
		}
	if (inValidCommand == YES)
		{
		for (i=0; i<MAX_NUM_DIVS; i++)
			activeParts[i] = NO;
		inValidCommand = NO;
		}

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		/* set Applyto (Applyto) *************************************************************/
		if (!strcmp(parmName, "Applyto"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(LEFTPAR);
			else if (expecting == Expecting(LEFTPAR))
				{
				for (i=0; i<MAX_NUM_DIVS; i++)
					activeParts[i] = NO;
				fromI = toJ = -1;
				foundDash = NO;
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				if (fromI != -1)
					activeParts[fromI-1] = YES;
#				if 0
				for (i=0; i<numCurrentDivisions; i++)
					MrBayesPrint("%d ", activeParts[i]);
				MrBayesPrint ("\n");
#				endif
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else if (expecting == Expecting(COMMA))
				{
				foundComma = YES;
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (IsSame ("All", tkn) == DIFFERENT)
					{
					MrBayesPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
					return (ERROR);
					}
				for (i=0; i<numCurrentDivisions; i++)
					activeParts[i] = YES;
				expecting  = Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt > numCurrentDivisions)
					{
					MrBayesPrint ("%s   Partition delimiter is too large\n", spacer);
					return (ERROR);
					}
				if (fromI == -1)
					fromI = tempInt;
				else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
					{
					toJ = tempInt;
					for (i=fromI-1; i<toJ; i++)
						activeParts[i] = YES;
					fromI = toJ = -1;
					foundDash = NO;
					}
				else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
					{
					activeParts[fromI-1] = YES;
					fromI = tempInt;
					foundComma = NO;
					}
					
				expecting  = Expecting(COMMA);
				expecting |= Expecting(DASH);
				expecting |= Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}
		/* set Tratiopr (tRatioPr) ************************************************************/
		else if (!strcmp(parmName, "Tratiopr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].tRatioPr, tempStr);
							modelParams[i].tRatioDir[0] = modelParams[i].tRatioDir[1] = 1.0;
							modelParams[i].tRatioFix = 1.0;
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Tratiopr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].tRatioPr,"Beta"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (tempD > ALPHA_MAX)
								{
								MrBayesPrint ("%s   Beta parameter cannot be greater than %1.2lf\n", spacer, ALPHA_MAX);
								return (ERROR);
								}
							if (tempD < ALPHA_MIN)
								{
								MrBayesPrint ("%s   Beta parameter cannot be less than %1.2lf\n", spacer, ALPHA_MIN);
								return (ERROR);
								}
							modelParams[i].tRatioDir[numVars[i]++] = tempD;
							if (numVars[i] < 2)
								expecting = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Tratiopr to Beta(%1.2lf,%1.2lf)\n", spacer, modelParams[i].tRatioDir[0], modelParams[i].tRatioDir[1]);
								else
									MrBayesPrint ("%s   Setting Tratiopr to Beta(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].tRatioDir[0], modelParams[i].tRatioDir[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].tRatioPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].tRatioFix = tempD;
							if (modelParams[i].tRatioFix > KAPPA_MAX)
								{
								MrBayesPrint ("%s   Tratio cannot be greater than %1.2lf\n", spacer, KAPPA_MAX);
								return (ERROR);
								}
							if (modelParams[i].tRatioFix < 0.0)
								{
								MrBayesPrint ("%s   Tratio cannot be less than %1.2lf\n", spacer, 0.0);
								return (ERROR);
								}
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Tratiopr to Fixed(%1.2lf)\n", spacer, modelParams[i].tRatioFix);
							else
								MrBayesPrint ("%s   Setting Tratiopr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].tRatioFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Revmatpr (revMatPr) ************************************************************/
		else if (!strcmp(parmName, "Revmatpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].revMatPr, tempStr);
							modelParams[i].revMatDir[0] = modelParams[i].revMatDir[1] = 1.0;
							modelParams[i].revMatDir[2] = modelParams[i].revMatDir[3] = 1.0;
							modelParams[i].revMatDir[4] = modelParams[i].revMatDir[5] = 1.0;
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Revmatpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				tempNumStates = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				/* find out what type of prior is being set */
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						strcpy (tempStr,modelParams[i].revMatPr);
					}
				/* find and store the number */
				sscanf (tkn, "%lf", &tempD);
				if (!strcmp(tempStr,"Dirichlet"))
					{
					if (tempD > ALPHA_MAX)
						{
						MrBayesPrint ("%s   Dirichlet parameter cannot be greater than %1.2lf\n", spacer, ALPHA_MAX);
						return (ERROR);
						}
					if (tempD < ALPHA_MIN)
						{
						MrBayesPrint ("%s   Dirichlet parameter cannot be less than %1.2lf\n", spacer, ALPHA_MIN);
						return (ERROR);
						}
					}
				else if (!strcmp(tempStr,"Fixed"))
					{
					if (tempD > KAPPA_MAX)
						{
						MrBayesPrint ("%s   Rate value cannot be greater than %1.2lf\n", spacer, KAPPA_MAX);
						return (ERROR);
						}
					if (tempD < 0.01)
						{
						MrBayesPrint ("%s   Rate value cannot be less than %1.2lf\n", spacer, 0.01);
						return (ERROR);
						}
					}
				tempNum[tempNumStates++] = tempD;
				if (tempNumStates == 1 && !strcmp(tempStr,"Dirichlet"))
					expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
				else if (tempNumStates < 6)
					expecting  = Expecting(COMMA);
				else
					expecting = Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].revMatPr,"Dirichlet"))
							{
							for (j=0; j<6; j++)
								{
								if (tempNumStates == 1)
									modelParams[i].revMatDir[j] = tempNum[0] / (MrBFlt) 6.0;
								else
									modelParams[i].revMatDir[j] = tempNum[j];
								}

							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Revmatpr to Dirichlet(%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer, 
								modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
								modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5]);
							else
								MrBayesPrint ("%s   Setting Revmatpr to Dirichlet(%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf) for partition %d\n", spacer, 
								modelParams[i].revMatDir[0], modelParams[i].revMatDir[1], modelParams[i].revMatDir[2],
								modelParams[i].revMatDir[3], modelParams[i].revMatDir[4], modelParams[i].revMatDir[5], i+1);
							}
						else if (!strcmp(modelParams[i].revMatPr,"Fixed"))
							{
							for (j=0; j<6; j++)
								modelParams[i].revMatFix[j] = tempNum[j];
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Revmatpr to Fixed(%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer, 
								modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2],
								modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5]);
							else
								MrBayesPrint ("%s   Setting Revmatpr to Fixed(%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf) for partition %d\n", spacer, 
								modelParams[i].revMatFix[0], modelParams[i].revMatFix[1], modelParams[i].revMatFix[2],
								modelParams[i].revMatFix[3], modelParams[i].revMatFix[4], modelParams[i].revMatFix[5], i+1);
							}
						}
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Aarevmatpr (aaRevMatPr) ********************************************************/
		else if (!strcmp(parmName, "Aarevmatpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == PROTEIN)
							{
							strcpy(modelParams[i].aaRevMatPr, tempStr);
							for (j=0; j<190; j++)
								modelParams[i].aaRevMatDir[j] = 1.0;
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Aarevmatpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				tempNumStates = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				/* find out what type of prior is being set */
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == PROTEIN)
						strcpy (tempStr,modelParams[i].aaRevMatPr);
					}
				/* find and store the number */
				sscanf (tkn, "%lf", &tempD);
				if (!strcmp(tempStr,"Dirichlet"))
					{
					if (tempD > ALPHA_MAX)
						{
						MrBayesPrint ("%s   Dirichlet parameter cannot be greater than %1.2lf\n", spacer, ALPHA_MAX);
						return (ERROR);
						}
					if (tempD < ALPHA_MIN)
						{
						MrBayesPrint ("%s   Dirichlet parameter cannot be less than %1.2lf\n", spacer, ALPHA_MIN);
						return (ERROR);
						}
					}
				else if (!strcmp(tempStr,"Fixed"))
					{
					if (tempD > KAPPA_MAX)
						{
						MrBayesPrint ("%s   Rate value cannot be greater than %1.2lf\n", spacer, KAPPA_MAX);
						return (ERROR);
						}
					if (tempD < 0.01)
						{
						MrBayesPrint ("%s   Rate value cannot be less than %1.2lf\n", spacer, 0.01);
						return (ERROR);
						}
					}
				tempStateFreqs[tempNumStates++] = tempD;
				if (tempNumStates == 1 && !strcmp(tempStr,"Dirichlet"))
					expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
				else if (tempNumStates < 190)
					expecting  = Expecting(COMMA);
				else
					expecting = Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == PROTEIN)
						{
						if (!strcmp(modelParams[i].aaRevMatPr,"Dirichlet"))
							{
							for (j=0; j<190; j++)
								{
								if (tempNumStates == 1)
									modelParams[i].aaRevMatDir[j] = tempStateFreqs[0] / (MrBFlt) 190.0;
								else
									modelParams[i].aaRevMatDir[j] = tempStateFreqs[j];
								}
							if (nApplied == 0 && numCurrentDivisions == 1)
								{
								for (j=0; j<190; j++)
									if (AreDoublesEqual(modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[j], 0.00001) == NO)
										break;
								if (j == 190)
									MrBayesPrint ("%s   Setting Aarevmatpr to Dirichlet(%1.2lf,%1.2lf,...)\n", spacer,
										modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[0]);
								else
									{
									MrBayesPrint ("%s   Setting Aarevmatpr to Dirichlet(\n", spacer);
									for (j=0; j<190; j++)
										{
										if (j % 10 == 0)
											MrBayesPrint ("%s      ", spacer);
										MrBayesPrint ("%1.2lf", modelParams[i].aaRevMatDir[j]);
										if (j == 189)
											MrBayesPrint (")\n");
										else if ((j+1) % 10 == 0)
											MrBayesPrint (",\n");
										else
											MrBayesPrint (",");
										}
									}
								}
							else
								{
								for (j=0; j<190; j++)
									if (AreDoublesEqual(modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[j], 0.00001) == NO)
										break;
								if (j == 190)
									MrBayesPrint ("%s   Setting Aarevmatpr to Dirichlet(%1.2lf,%1.2lf,...) for partition %d\n",
										spacer, modelParams[i].aaRevMatDir[0], modelParams[i].aaRevMatDir[0], i+1);
								else
									{
									MrBayesPrint ("%s   Setting Aarevmatpr to Dirichlet(\n", spacer);
									for (j=0; j<190; j++)
										{
										if (j % 10 == 0)
											MrBayesPrint ("%s      ", spacer);
										MrBayesPrint ("%1.2lf", modelParams[i].aaRevMatDir[j]);
										if (j == 189)
											MrBayesPrint (")\n");
										else if ((j+1) % 10 == 0)
											MrBayesPrint (",\n");
										else
											MrBayesPrint (",");
										}
									}
									MrBayesPrint ("%s      for partition %d\n", spacer, i+1);
								}
							}
						else if (!strcmp(modelParams[i].aaRevMatPr,"Fixed"))
							{
							for (j=0; j<190; j++)
								modelParams[i].aaRevMatFix[j] = tempStateFreqs[j];
							if (nApplied == 0 && numCurrentDivisions == 1)
								{
								for (j=0; j<190; j++)
									if (AreDoublesEqual(modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[j], 0.00001) == NO)
										break;
								if (j == 190)
									MrBayesPrint ("%s   Setting Aarevmatpr to Fixed(%1.2lf,%1.2lf,...)\n", spacer, modelParams[i].aaRevMatFix[0],
										modelParams[i].aaRevMatFix[0]);
								else
									{
									MrBayesPrint ("%s   Setting Aarevmatpr to Fixed(\n", spacer);
									for (j=0; j<190; j++)
										{
										if (j % 10 == 0)
											MrBayesPrint ("%s      ", spacer);
										MrBayesPrint ("%1.2lf", modelParams[i].aaRevMatFix[j]);
										if (j == 189)
											MrBayesPrint (")\n");
										else if ((j+1) % 10 == 0)
											MrBayesPrint (",\n");
										else
											MrBayesPrint (",");
										}
									}
								}
							else
								{
								for (j=0; j<190; j++)
									if (AreDoublesEqual(modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[j], 0.00001) == NO)
										break;
								if (j == 190)
									MrBayesPrint ("%s   Setting Aarevmatpr to Fixed(%1.2lf,%1.2lf,...) for partition %d\n", spacer,
										modelParams[i].aaRevMatFix[0], modelParams[i].aaRevMatFix[0], i+1);
								else
									{
									MrBayesPrint ("%s   Setting Aarevmatpr to Fixed(\n", spacer);
									for (j=0; j<190; j++)
										{
										if (j % 10 == 0)
											MrBayesPrint ("%s      ", spacer);
										MrBayesPrint ("%1.2lf", modelParams[i].aaRevMatFix[j]);
										if (j == 189)
											MrBayesPrint (")\n");
										else if ((j+1) % 10 == 0)
											MrBayesPrint (",\n");
										else
											MrBayesPrint (",");
										}
									}
									MrBayesPrint ("%s      for partition %d\n", spacer, i+1);
								}
							}
						}
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Omegapr (omegaPr) **************************************************************/
		else if (!strcmp(parmName, "Omegapr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].omegaPr, tempStr);
							modelParams[i].omegaDir[0] = modelParams[i].omegaDir[1] = 1.0;
							modelParams[i].omegaFix = 1.0;
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Omegapr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].omegaPr,"Dirichlet"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (tempD > ALPHA_MAX)
								{
								MrBayesPrint ("%s   Dirichlet parameter cannot be greater than %1.2lf\n", spacer, ALPHA_MAX);
								return (ERROR);
								}
							if (tempD < ALPHA_MIN)
								{
								MrBayesPrint ("%s   Dirichlet parameter cannot be less than %1.2lf\n", spacer, ALPHA_MIN);
								return (ERROR);
								}
							modelParams[i].omegaDir[numVars[i]++] = tempD;
							if (numVars[i] < 1)
								expecting = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Omegapr to Dirichlet(%1.2lf,%1.2lf)\n", spacer, modelParams[i].omegaDir[0], modelParams[i].omegaDir[1]);
								else
									MrBayesPrint ("%s   Setting Omegapr to Dirichlet(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].omegaDir[0], modelParams[i].omegaDir[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].omegaPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].omegaFix = tempD;
							if (modelParams[i].omegaFix > KAPPA_MAX)
								{
								MrBayesPrint ("%s   Omega ratio cannot be greater than %1.2lf\n", spacer, KAPPA_MAX);
								return (ERROR);
								}
							if (modelParams[i].omegaFix < 0.0)
								{
								MrBayesPrint ("%s   Omega ratio cannot be less than %1.2lf\n", spacer, 0.0);
								return (ERROR);
								}
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Omegapr to Fixed(%1.2lf)\n", spacer, modelParams[i].omegaFix);
							else
								MrBayesPrint ("%s   Setting Omegapr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].omegaFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Ny98omega1pr (ny98omega1pr) ********************************************************/
		else if (!strcmp(parmName, "Ny98omega1pr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].ny98omega1pr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Ny98omega1pr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].ny98omega1pr,"Beta"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].ny98omega1Beta[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].ny98omega1Beta[0] < 0 || modelParams[i].ny98omega1Beta[1] < 0)
									{
									MrBayesPrint ("%s   Beta parameter should be greater than 0\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Ny98omega1pr to Beta(%1.2lf,%1.2lf)\n", spacer, modelParams[i].ny98omega1Beta[0], modelParams[i].ny98omega1Beta[1]);
								else
									MrBayesPrint ("%s   Setting Ny98omega1pr to Beta(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].ny98omega1Beta[0], modelParams[i].ny98omega1Beta[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].ny98omega1pr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].ny98omega1Fixed = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Ny98omega1pr to Fixed(%1.2lf)\n", spacer, modelParams[i].ny98omega1Fixed);
							else
								MrBayesPrint ("%s   Setting Ny98omega1pr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].ny98omega1Fixed, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Ny98omega3pr (ny98omega3pr) ********************************************************/
		else if (!strcmp(parmName, "Ny98omega3pr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].ny98omega3pr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Ny98omega3pr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].ny98omega3pr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].ny98omega3Uni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].ny98omega3Uni[0] >= modelParams[i].ny98omega3Uni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Ny98omega3pr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].ny98omega3Uni[0], modelParams[i].ny98omega3Uni[1]);
								else
									MrBayesPrint ("%s   Setting Ny98omega3pr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].ny98omega3Uni[0], modelParams[i].ny98omega3Uni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].ny98omega3pr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].ny98omega3Exp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Ny98omega3pr to Exponential(%1.2lf)\n", spacer, modelParams[i].ny98omega3Exp);
							else
								MrBayesPrint ("%s   Setting Ny98omega3pr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].ny98omega3Exp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].ny98omega3pr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].ny98omega3Fixed = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Ny98omega3pr to Fixed(%1.2lf)\n", spacer, modelParams[i].ny98omega3Fixed);
							else
								MrBayesPrint ("%s   Setting Ny98omega3pr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].ny98omega3Fixed, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set M3omegapr (m3omegapr) ********************************************************/
		else if (!strcmp(parmName, "M3omegapr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].m3omegapr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid M3omegapr argument\n", spacer);
					return (ERROR);
					}
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].m3omegapr,"Exponential"))
							{
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting M3omegapr to Exponential\n", spacer);
							else
								MrBayesPrint ("%s   Setting M3omegapr to Exponential for partition %d\n", spacer, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				if (!strcmp(tempStr,"Exponential"))
					expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				else
					expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].m3omegapr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m3omegaFixed[numVars[i]++] = tempD;
							if (numVars[i] == 1 || numVars[i] == 2)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].m3omegaFixed[0] >= modelParams[i].m3omegaFixed[1] || modelParams[i].m3omegaFixed[0] >= modelParams[i].m3omegaFixed[2] || modelParams[i].m3omegaFixed[1] >= modelParams[i].m3omegaFixed[2])
									{
									MrBayesPrint ("%s   The three omega values must be ordered, such that omega1 < omega2 < omega3\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting M3omegapr to Fixed(%1.2lf,%1.2lf,%1.2lf)\n", spacer, modelParams[i].m3omegaFixed[0], modelParams[i].m3omegaFixed[1], modelParams[i].m3omegaFixed[2]);
								else
									MrBayesPrint ("%s   Setting M3omegapr to Fixed(%1.2lf,%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].m3omegaFixed[0], modelParams[i].m3omegaFixed[1], modelParams[i].m3omegaFixed[2], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Codoncatfreqs (codonCatFreqPr) ********************************************************/
		else if (!strcmp(parmName, "Codoncatfreqs"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].codonCatFreqPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Omegapurpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].codonCatFreqPr,"Dirichlet"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].codonCatDir[numVars[i]++] = tempD;
							if (numVars[i] == 1 || numVars[i] == 2)
								expecting  = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Codoncatfreqs prior to Dirichlet(%1.2lf,%1.2lf,%1.2lf)\n", spacer, modelParams[i].codonCatDir[0], modelParams[i].codonCatDir[1], modelParams[i].codonCatDir[2]);
								else
									MrBayesPrint ("%s   Setting Codoncatfreqs prior to Dirichlet(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].codonCatDir[0], modelParams[i].codonCatDir[1], modelParams[i].codonCatDir[2], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].codonCatFreqPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].codonCatFreqFix[numVars[i]++] = tempD;
							if (numVars[i] == 1 || numVars[i] == 2)
								expecting  = Expecting(COMMA);
							else
								{
								if (AreDoublesEqual (modelParams[i].codonCatFreqFix[0] + modelParams[i].codonCatFreqFix[1] + modelParams[i].codonCatFreqFix[2], (MrBFlt) 1.0, (MrBFlt) 0.001) == NO)
									{
									MrBayesPrint ("%s   Codon category frequencies must sum to 1\n", spacer);
									return (ERROR);
									}
								
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Codoncatfreqs prior to Fixed(%1.2lf,%1.2lf,%1.2lf)\n", spacer, modelParams[i].codonCatFreqFix[0], modelParams[i].codonCatFreqFix[1], modelParams[i].codonCatFreqFix[2]);
								else
									MrBayesPrint ("%s   Setting Codoncatfreqs prior to Fixed(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].codonCatFreqFix[0], modelParams[i].codonCatFreqFix[1], modelParams[i].codonCatFreqFix[2], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}

		/* set Shapepr (shapePr) **************************************************************/
		else if (!strcmp(parmName, "Shapepr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN || modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD))
							strcpy(modelParams[i].shapePr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Shapepr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN || modelParams[i].dataType == RESTRICTION || modelParams[i].dataType == STANDARD))
						{
						if (!strcmp(modelParams[i].shapePr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].shapeUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].shapeUni[0] >= modelParams[i].shapeUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (modelParams[i].shapeUni[1] > MAX_SHAPE_PARAM)
									{
									MrBayesPrint ("%s   Upper value for uniform cannot be greater than %1.2lf\n", spacer, MAX_SHAPE_PARAM);
									return (ERROR);
									}
								if (modelParams[i].shapeUni[0] < MIN_SHAPE_PARAM)
									{
									MrBayesPrint ("%s   Lower value for uniform cannot be less than %1.2lf\n", spacer, MIN_SHAPE_PARAM);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Shapepr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1]);
								else
									MrBayesPrint ("%s   Setting Shapepr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].shapeUni[0], modelParams[i].shapeUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].shapePr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].shapeExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Shapepr to Exponential(%1.2lf)\n", spacer, modelParams[i].shapeExp);
							else
								MrBayesPrint ("%s   Setting Shapepr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].shapeExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].shapePr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].shapeFix = tempD;
							if (modelParams[i].shapeFix > MAX_SHAPE_PARAM)
								{
								MrBayesPrint ("%s   Shape parameter cannot be greater than %1.2lf\n", spacer, MAX_SHAPE_PARAM);
								return (ERROR);
								}
							if (modelParams[i].shapeFix < MIN_SHAPE_PARAM)
								{
								MrBayesPrint ("%s   Shape parameter cannot be less than %1.2lf\n", spacer, MIN_SHAPE_PARAM);
								return (ERROR);
								}
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Shapepr to Fixed(%1.2lf)\n", spacer, modelParams[i].shapeFix);
							else
								MrBayesPrint ("%s   Setting Shapepr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].shapeFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Pinvarpr (pInvarPr) ************************************************************/
		else if (!strcmp(parmName, "Pinvarpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].pInvarPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Pinvarpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].pInvarPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].pInvarUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].pInvarUni[0] >= modelParams[i].pInvarUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (modelParams[i].pInvarUni[1] > 1.0)
									{
									MrBayesPrint ("%s   Upper value for uniform should be less than or equal to 1.0\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Pinvarpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1]);
								else
									MrBayesPrint ("%s   Setting Pinvarpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].pInvarUni[0], modelParams[i].pInvarUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].pInvarPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (tempD > 1.0)
								{
								MrBayesPrint ("%s   Value for Pinvar should be in the interval (0, 1)\n", spacer);
								return (ERROR);
								}
							modelParams[i].pInvarFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Pinvarpr to Fixed(%1.2lf)\n", spacer, modelParams[i].pInvarFix);
							else
								MrBayesPrint ("%s   Setting Pinvarpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].pInvarFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Ratecorrpr (adGammaCorPr) ******************************************************/
		else if (!strcmp(parmName, "Ratecorrpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].adGammaCorPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Ratecorrpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				foundDash = NO;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting = Expecting(NUMBER) | Expecting(DASH);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].adGammaCorPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (foundDash == YES)
							tempD *= -1.0;
							modelParams[i].corrUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].corrUni[0] >= modelParams[i].corrUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (modelParams[i].corrUni[1] > 1.0)
									{
									MrBayesPrint ("%s   Upper value for uniform should be less than or equal to 1.0\n", spacer);
									return (ERROR);
									}
								if (modelParams[i].corrUni[0] < -1.0)
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than or equal to -1.0\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Ratecorrpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].corrUni[0], modelParams[i].corrUni[1]);
								else
									MrBayesPrint ("%s   Setting Ratecorrpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].corrUni[0], modelParams[i].corrUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].adGammaCorPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (foundDash == YES)
								tempD *= -1.0;
							if (tempD > 1.0 || tempD < -1.0)
								{
								MrBayesPrint ("%s   Value for Ratecorrpr should be in the interval (-1, +1)\n", spacer);
								return (ERROR);
								}
							modelParams[i].corrFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Ratecorrpr to Fixed(%1.2lf)\n", spacer, modelParams[i].corrFix);
							else
								MrBayesPrint ("%s   Setting Ratecorrpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].corrFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				foundDash = NO;
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER) | Expecting(DASH);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Browncorrpr (brownCorPr) ******************************************************/
		else if (!strcmp(parmName, "Browncorrpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
							strcpy(modelParams[i].brownCorPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Browncorrpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				foundDash = NO;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting = Expecting(NUMBER) | Expecting(DASH);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
						{
						if (!strcmp(modelParams[i].brownCorPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (foundDash == YES)
							tempD *= -1.0;
							modelParams[i].brownCorrUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].brownCorrUni[0] >= modelParams[i].brownCorrUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (modelParams[i].brownCorrUni[1] > 1.0)
									{
									MrBayesPrint ("%s   Upper value for uniform should be less than or equal to 1.0\n", spacer);
									return (ERROR);
									}
								if (modelParams[i].brownCorrUni[0] < -1.0)
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than or equal to -1.0\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Browncorrpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownCorrUni[0], modelParams[i].brownCorrUni[1]);
								else
									MrBayesPrint ("%s   Setting Browncorrpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brownCorrUni[0], modelParams[i].brownCorrUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].brownCorPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (foundDash == YES)
								tempD *= -1.0;
							if (tempD > 1.0 || tempD < -1.0)
								{
								MrBayesPrint ("%s   Value for Browncorrpr should be in the interval (-1, +1)\n", spacer);
								return (ERROR);
								}
							modelParams[i].brownCorrFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Browncorrpr to Fixed(%1.2lf)\n", spacer, modelParams[i].brownCorrFix);
							else
								MrBayesPrint ("%s   Setting Browncorrpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].brownCorrFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				foundDash = NO;
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER) | Expecting(DASH);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Ratepr (ratePr) *****************************************************************/
		else if (!strcmp(parmName, "Ratepr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
							{
							if (!strcmp(tempStr,"Variable"))
								strcpy(modelParams[i].ratePr, "Dirichlet");
							else
								strcpy(modelParams[i].ratePr, tempStr);
							modelParams[i].ratePrDir = 1.0;
							if (!strcmp(tempStr,"Variable") || !strcmp(tempStr,"Fixed"))
								{
								if (tempStr[0]=='V')
									strcat (tempStr," [Dirichlet(..,1,..)]");
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Ratepr to %s\n", spacer, tempStr);
								else
									MrBayesPrint ("%s   Setting Ratepr to %s for partition %d\n", spacer, tempStr, i+1);
								if (tempStr[0]=='V')
									strcpy (tempStr,"Variable");
								}
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Ratepr argument\n", spacer);
					return (ERROR);
					}
				if (!strcmp(tempStr,"Fixed") || !strcmp(tempStr,"Variable"))
					expecting  = Expecting(PARAMETER) | Expecting(SEMICOLON);
				else
					expecting = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting = Expecting (NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				/* find next partition to fill in */
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					if ((activeParts[i] == YES || nApplied == 0) && numVars[i] == 0)
						break;
				if (i == numCurrentDivisions)
					{
					MrBayesPrint ("%s   Could not find first ratemultiplier partition\n", spacer);
					return (ERROR);
					}
				numVars[i] = 1;
				/* read in the parameter */
				sscanf (tkn, "%lf", &tempD);
				if (tempD < ALPHA_MIN || tempD > ALPHA_MAX)
					{
					MrBayesPrint ("%s   Ratemultiplier Dirichlet parameter %lf out of range\n", spacer, tempD);
					return (ERROR);
					}
				/* set the parameter */
				modelParams[i].ratePrDir = tempD;				
				/* check if all partitions have been filled in */
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && numVars[i] == 0)
						break;
					}
				/* set expecting accordingly so that we know what should be coming next */
				if (i == numCurrentDivisions)
					expecting = Expecting (RIGHTPAR);
				else
					expecting = Expecting (COMMA);
				}
			else if (expecting == Expecting (COMMA))
				expecting = Expecting (NUMBER);
			else if (expecting == Expecting (RIGHTPAR))
				{
				/* print message */
				for (i=j=0; i<numCurrentDivisions; i++)
					{
					if (numVars[i] == 1)
						{
						j++;
						if (j == 1)
							{
							MrBayesPrint ("%s   Setting Ratepr to Dirichlet(%1.2f",
								spacer, modelParams[i].ratePrDir);
							}
						else
							MrBayesPrint(",%1.2f", modelParams[i].ratePrDir);
						}
					}
				if (numCurrentDivisions == 1)
					MrBayesPrint (")\n");
				else
					{
					MrBayesPrint (") for partition");
					if (j > 1)
						MrBayesPrint ("s");
					for (i=k=0; i<numCurrentDivisions; i++)
						{
						if (numVars[i] == 1)
							{
							k++;
							if (k == j && j > 1)
								MrBayesPrint (", and %d", i+1);
							else if (k == 1)
								MrBayesPrint (" %d", i+1);
							else
								MrBayesPrint (", %d", i+1);
							}
						}
					MrBayesPrint ("\n");
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Covswitchpr (covSwitchPr) ******************************************************/
		else if (!strcmp(parmName, "Covswitchpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
							strcpy(modelParams[i].covSwitchPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Covswitchpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
						{
						if (!strcmp(modelParams[i].covSwitchPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].covswitchUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].covswitchUni[0] >= modelParams[i].covswitchUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Covswitchpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1]);
								else
									MrBayesPrint ("%s   Setting Covswitchpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].covswitchUni[0], modelParams[i].covswitchUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].covSwitchPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].covswitchExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Covswitchpr to Exponential(%1.2lf)\n", spacer, modelParams[i].covswitchExp);
							else
								MrBayesPrint ("%s   Setting Covswitchpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].covswitchExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].covSwitchPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].covswitchFix[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Covswitchpr to Fixed(%1.4lf,%1.4lf)\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[1]);
								else
									MrBayesPrint ("%s   Setting Covswitchpr to Fixed(%1.4lf,%1.4lf) for partition %d\n", spacer, modelParams[i].covswitchFix[0], modelParams[i].covswitchFix[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Symdirihyperpr (symPiPr) ******************************************************/
		else if (!strcmp(parmName, "Symdirihyperpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				foundBeta = NO;
				expecting = Expecting(ALPHA);
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (foundBeta == NO)
					{
					/* expecting to see Uniform, Exponential, or Fixed */
					if (IsArgValid(tkn, tempStr) == NO_ERROR)
						{
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
								{
								strcpy(modelParams[i].symPiPr, tempStr);
								}
							}
						}
					else
						{
						MrBayesPrint ("%s   Invalid Symdirihyperpr argument\n", spacer);
						return (ERROR);
						}
					expecting  = Expecting(LEFTPAR);
					for (i=0; i<MAX_NUM_DIVS; i++)
						numVars[i] = 0;
					foundBeta = YES;	
					}	
				else
					{
					/* expecting infinity */
					if (IsSame("Infinity", tkn) == SAME || IsSame("Infinity", tkn) == CONSISTENT_WITH)
						{
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
								{
								if (!strcmp(modelParams[i].symPiPr, "Fixed"))
									{
									modelParams[i].symBetaFix = -1;
									if (nApplied == 0 && numCurrentDivisions == 1)
										MrBayesPrint ("%s   Setting Symdirihyperpr to Beta(Infinity)\n", spacer);
									else
										MrBayesPrint ("%s   Setting Symdirihyperpr to Beta(Infinity) for partition %d\n", spacer, i+1);
									expecting  = Expecting(RIGHTPAR);
									}
								else
									{
									MrBayesPrint ("%s   Problem setting Symdirihyperpr\n", spacer);
									return (ERROR);
									}
								}
							}
						expecting  = Expecting(RIGHTPAR);
						}
					}		
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == STANDARD || modelParams[i].dataType == RESTRICTION))
						{
						if (!strcmp(modelParams[i].symPiPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].symBetaFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Symdirihyperpr to Fixed(%1.2lf)\n", spacer, modelParams[i].symBetaFix);
							else
								MrBayesPrint ("%s   Setting Symdirihyperpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].symBetaFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].symPiPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].symBetaExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Symdirihyperpr to Exponential(%1.2lf)\n", spacer, modelParams[i].symBetaExp);
							else
								MrBayesPrint ("%s   Setting Symdirihyperpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].symBetaExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].symPiPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].symBetaUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)	
								expecting = Expecting(COMMA);
							else
								{
								if (modelParams[i].symBetaUni[0] >= modelParams[i].symBetaUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Symdirihyperpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1]);
								else
									MrBayesPrint ("%s   Setting Symdirihyperpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].symBetaUni[0], modelParams[i].symBetaUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else
							{
							MrBayesPrint ("%s   Problem setting Symdirihyperpr\n", spacer);
							return (ERROR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Statefreqpr (stateFreqPr) ******************************************************/
		else if (!strcmp(parmName, "Statefreqpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsSame ("Equal", tkn) == DIFFERENT && IsSame ("Empirical", tkn) == DIFFERENT)
					{
					/* the user wants to specify a dirichlet or fixed prior */
					if (IsArgValid(tkn, tempStr) == NO_ERROR)
						{
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
								strcpy(modelParams[i].stateFreqPr, tempStr);
							}
						}
					else
						{
						MrBayesPrint ("%s   Invalid Statefreqpr argument\n", spacer);
						return (ERROR);
						}
					expecting  = Expecting(LEFTPAR);
					}
				else
					{
					/* the user wants equal or empirical state frequencies */
					nApplied = NumActiveParts ();
					if (IsSame ("Equal", tkn) == SAME || IsSame ("Equal", tkn) == CONSISTENT_WITH)
						{
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
								strcpy(modelParams[i].stateFreqsFixType, "Equal");
							}
						}
					else if (IsSame ("Empirical", tkn) == SAME || IsSame ("Empirical", tkn) == CONSISTENT_WITH)
						{
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
								strcpy(modelParams[i].stateFreqsFixType, "Empirical");
							}
						}
					else
						{
						MrBayesPrint ("%s   Invalid Statefreqpr delimiter\n", spacer);
						return (ERROR);
						}
					expecting  = Expecting(RIGHTPAR);
					}
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				tempNumStates = 0;
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				sscanf (tkn, "%lf", &tempD);
				tempStateFreqs[tempNumStates++] = tempD;
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
						if (!strcmp(modelParams[i].stateFreqPr,"Fixed"))
							strcpy(modelParams[i].stateFreqsFixType, "User");
					}
				expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType != CONTINUOUS)
						{
						ns = NumStates(i);
						if (!strcmp(modelParams[i].stateFreqPr,"Dirichlet"))
							{
							if (tempNumStates == 1)
								{
								for (j=0; j<ns; j++)
									modelParams[i].stateFreqsDir[j] = tempStateFreqs[0] / ns;
								MrBayesPrint ("%s   Setting Statefreqpr to Dirichlet(", spacer);
								for (j=0; j<ns; j++)
									{
									MrBayesPrint("%1.2lf", modelParams[i].stateFreqsDir[j]);
									if (j == ns - 1)
										MrBayesPrint (")");
									else
										MrBayesPrint (",");	
									}	
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("\n");
								else
									MrBayesPrint (" for partition %d\n", i+1); 
								modelParams[i].numDirParams = ns;
								}
							else
								{
								if (tempNumStates != ns)
									{
									MrBayesPrint ("%s   Found %d dirichlet parameters but expecting %d\n", spacer, tempNumStates, modelParams[i].nStates);
									return (ERROR);
									}
								else
									{
									modelParams[i].numDirParams = ns;
									for (j=0; j<ns; j++)
										modelParams[i].stateFreqsDir[j] = tempStateFreqs[j];
									MrBayesPrint ("%s   Setting Statefreqpr to Dirichlet(", spacer);
									for (j=0; j<ns; j++)
										{
										MrBayesPrint("%1.2lf", modelParams[i].stateFreqsDir[j]);
										if (j == ns - 1)
											MrBayesPrint (")");
										else
											MrBayesPrint (",");	
										}	
									if (nApplied == 0 && numCurrentDivisions == 1)
										MrBayesPrint ("\n");
									else
										MrBayesPrint (" for partition %d\n", i+1); 
									}
								}
							}
						else if (!strcmp(modelParams[i].stateFreqPr,"Fixed"))
							{
							if (tempNumStates == 0)
								{
								if (!strcmp(modelParams[i].stateFreqsFixType, "Equal"))
									MrBayesPrint ("%s   Setting Statefreqpr to Fixed(Equal)", spacer);
								else if (!strcmp(modelParams[i].stateFreqsFixType, "Empirical"))
									MrBayesPrint ("%s   Setting Statefreqpr to Fixed(Empirical)", spacer);
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("\n");
								else
									MrBayesPrint (" for partition %d\n", i+1); 
								}
							else 
								{
								if (tempNumStates == ns)
									{
									sum = 0.0;
									for (j=0; j<ns; j++)
										sum += tempStateFreqs[j];
									if (AreDoublesEqual (sum, (MrBFlt) 1.0, (MrBFlt) 0.001) == NO)
										{
										MrBayesPrint ("%s   State frequencies do not sum to 1.0\n", spacer);
										return (ERROR);
										}
									strcpy(modelParams[i].stateFreqsFixType, "User");
									for (j=0; j<ns; j++)
										modelParams[i].stateFreqsFix[j] = tempStateFreqs[j];
									MrBayesPrint ("%s   Setting Statefreqpr to Fixed(", spacer);
									for (j=0; j<ns; j++)
										{
										MrBayesPrint("%1.2lf", modelParams[i].stateFreqsFix[j]);
										if (j == ns - 1)
											MrBayesPrint (")");
										else
											MrBayesPrint (",");	
										}	
									if (nApplied == 0 && numCurrentDivisions == 1)
										MrBayesPrint ("\n");
									else
										MrBayesPrint (" for partition %d\n", i+1); 
									}
								else
									{
									MrBayesPrint ("%s   Found %d state frequencies but expecting %d\n", spacer, tempNumStates, modelParams[i].nStates);
									return (ERROR);
									}
								}
								
							}
						}
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Topologypr (topologyPr) ********************************************************/
		else if (!strcmp(parmName, "Topologypr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				foundEqual = YES;
				expecting = Expecting(ALPHA);
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (foundEqual == YES)
					{
					if (IsArgValid(tkn, tempStr) == NO_ERROR)
						{
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if (activeParts[i] == YES || nApplied == 0)
								{
								strcpy(modelParams[i].topologyPr, tempStr);
								if (!strcmp(modelParams[i].topologyPr, "Constraints"))
									for (j=0; j<MAX_NUM_CONSTRAINTS; j++)
										modelParams[i].activeConstraints[j] = NO;
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Topologypr to %s\n", spacer, modelParams[i].topologyPr);
								else
									MrBayesPrint ("%s   Setting Topologypr to %s for partition %d\n", spacer, modelParams[i].topologyPr, i+1);
								}
							}
						}
					else
						{
						MrBayesPrint ("%s   Invalid Topologypr argument\n", spacer);
						return (ERROR);
						}
					if (!strcmp(tempStr, "Uniform"))
						expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
					else
						expecting = Expecting(LEFTPAR);
					foundEqual = NO;
					}
				else
					{
					if (foundDash == YES)
						{
						MrBayesPrint ("%s   Expecting a number\n", spacer);
						return (ERROR);
						}
					if (CheckString (tkn, constraintNames, &howMany) == ERROR)
						{
						MrBayesPrint ("%s   Could not find constraint named %s\n", spacer, tkn);
						return (ERROR);
						}
					numVars[howMany - 1] = YES;
					expecting  = Expecting(RIGHTPAR);
					expecting |= Expecting(COMMA);
					}
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = NO;
				fromI = toJ = -1;
				foundDash = foundComma = NO;
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(NUMBER))
				{
				if (numDefinedConstraints == 0)
					{
					MrBayesPrint ("%s   No constraints have been defined\n", spacer);
					return (ERROR);
					}
				sscanf (tkn, "%d", &tempInt);
				if (tempInt > numDefinedConstraints)
					{
					MrBayesPrint ("%s   Constraint number is too large\n", spacer);
					return (ERROR);
					}
				if (fromI == -1)
					fromI = tempInt;
				else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
					{
					toJ = tempInt;
					for (i=fromI-1; i<toJ; i++)
						numVars[i] = YES;
					fromI = toJ = -1;
					foundDash = NO;
					}
				else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
					{
					numVars[fromI-1] = YES;
					fromI = tempInt;
					foundComma = NO;
					}
				expecting  = Expecting(COMMA);
				expecting |= Expecting(DASH);
				expecting |= Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(COMMA))
				{
				foundComma = YES;
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				if (fromI != -1)
					numVars[fromI-1] = YES;
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if (activeParts[i] == YES || nApplied == 0)
						{
						modelParams[i].numActiveConstraints = 0;
						for (j=0; j<MAX_NUM_CONSTRAINTS; j++)
							{
							if (numVars[j] == YES)
								{
								modelParams[i].activeConstraints[j] = YES;
								modelParams[i].numActiveConstraints++;
								}
							else
								modelParams[i].activeConstraints[j] = NO;
							}
						if (modelParams[i].numActiveConstraints == 0)
							{
							MrBayesPrint ("%s   No constraints have been defined\n", spacer);
							return (ERROR);
							}

						}
					}
#				if 0
				for (i=0; i<numCurrentDivisions; i++)
					{
					MrBayesPrint ("%4d -- ", i+1);
					for (j=0; j<numDefinedConstraints; j++)
						MrBayesPrint (" %d", modelParams[i].activeConstraints[j]);
					MrBayesPrint ("\n");
					}
#				endif				
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Brlenspr (brlensPr) ************************************************************/
		else if (!strcmp(parmName, "Brlenspr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = NO;
				foundEqual = YES;
				expecting = Expecting(ALPHA);
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (foundEqual == YES)
					{
					if (IsArgValid(tkn, tempStr) == NO_ERROR)
						{
						strcpy (colonPr, tempStr);
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							if (activeParts[i] == YES || nApplied == 0)
								strcpy(modelParams[i].brlensPr, tempStr);
						}
					else
						{
						MrBayesPrint ("%s   Invalid Brlenspr argument\n", spacer);
						return (ERROR);
						}
					foundEqual = NO;
					expecting  = Expecting(COLON);
					}
				else
					{
					if (!strcmp(colonPr, "Unconstrained"))
						{
						/* have unconstrained branch lengths, which we expect to have a uniform or exponential distribution */
						nApplied = NumActiveParts ();
						if (IsSame ("Uniform", tkn) == SAME || IsSame ("Uniform", tkn) == CONSISTENT_WITH)
							{
							for (i=0; i<numCurrentDivisions; i++)
								if (activeParts[i] == YES || nApplied == 0)
									strcpy(modelParams[i].unconstrainedPr, "Uniform");
							}
						else if (IsSame ("Exponential", tkn) == SAME || IsSame ("Exponential", tkn) == CONSISTENT_WITH)
							{
							for (i=0; i<numCurrentDivisions; i++)
								if (activeParts[i] == YES || nApplied == 0)
									strcpy(modelParams[i].unconstrainedPr, "Exponential");
							}
						else
							{
							MrBayesPrint ("%s   Do not understand %s\n", spacer, tkn);
							return (ERROR);
							}
						expecting  = Expecting(LEFTPAR);
						}
					else
						{
						/* otherwise we have a clock constraint and expect uniform, birthdeath, or coalescence prior */
						nApplied = NumActiveParts ();
						if (IsSame ("Uniform", tkn) == SAME || IsSame ("Uniform", tkn) == CONSISTENT_WITH)
							{
							for (i=0; i<numCurrentDivisions; i++)
								{
								if (activeParts[i] == YES || nApplied == 0)
									strcpy(modelParams[i].clockPr, "Uniform");
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Brlenspr to Clock:Uniform\n", spacer);
								else
									MrBayesPrint ("%s   Setting Brlenspr to Clock:Uniform for partition %d\n", spacer, i+1);
								}
							}
						else if (IsSame ("Birthdeath", tkn) == SAME || IsSame ("Birthdeath", tkn) == CONSISTENT_WITH)
							{
							for (i=0; i<numCurrentDivisions; i++)
								{
								if (activeParts[i] == YES || nApplied == 0)
									strcpy(modelParams[i].clockPr, "Birthdeath");
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Brlenspr to Clock:Birthdeath\n", spacer);
								else
									MrBayesPrint ("%s   Setting Brlenspr to Clock:Birthdeath for partition %d\n", spacer, i+1);
								}
							}
						else if (IsSame ("Coalescence", tkn) == SAME || IsSame ("Coalescence", tkn) == CONSISTENT_WITH)
							{
							for (i=0; i<numCurrentDivisions; i++)
								{
								if (activeParts[i] == YES || nApplied == 0)
									strcpy(modelParams[i].clockPr, "Coalescence");
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Brlenspr to Clock:Coalescence\n", spacer);
								else
									MrBayesPrint ("%s   Setting Brlenspr to Clock:Coalescence for partition %d\n", spacer, i+1);
								}
							}
						else
							{
							MrBayesPrint ("%s   Do not understand %s\n", spacer, tkn);
							return (ERROR);
							}
						expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
						}
					}
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				if (!strcmp(colonPr, "Unconstrained"))
					{
					/* have unconstrained branch lengths */
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							{
							if (!strcmp(modelParams[i].unconstrainedPr,"Uniform"))
								{
								sscanf (tkn, "%lf", &tempD);
								modelParams[i].brlensUni[numVars[i]++] = tempD;
								if (numVars[i] == 1)
									expecting  = Expecting(COMMA);
								else
									{
									if (modelParams[i].brlensUni[0] > 0.000001)
										{
										MrBayesPrint ("%s   Lower value for uniform must equal 0.0\n", spacer);
										return (ERROR);
										}
									if (modelParams[i].brlensUni[0] >= modelParams[i].brlensUni[1])
										{
										MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
										return (ERROR);
										}
									if (nApplied == 0 && numCurrentDivisions == 1)
										MrBayesPrint ("%s   Setting Brlenspr to Unconstrained:Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brlensUni[0], modelParams[i].brlensUni[1]);
									else
										MrBayesPrint ("%s   Setting Brlenspr to Unconstrained:Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brlensUni[0], modelParams[i].brlensUni[1], i+1);
									expecting  = Expecting(RIGHTPAR);
									}
								}
							else if (!strcmp(modelParams[i].unconstrainedPr,"Exponential"))
								{
								sscanf (tkn, "%lf", &tempD);
								modelParams[i].brlensExp = tempD;
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Brlenspr to Unconstrained:Exponential(%1.2lf)\n", spacer, modelParams[i].brlensExp);
								else
									MrBayesPrint ("%s   Setting Brlenspr to Unconstrained:Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].brlensExp, i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						}
					}
				else
					{
					
					}
				}
			else if (expecting == Expecting(COLON))
				{
				expecting  = Expecting(ALPHA);
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Speciationpr (speciationPr) ****************************************************/
		else if (!strcmp(parmName, "Speciationpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].speciationPr, tempStr);
					}
				else
					{
					MrBayesPrint ("%s   Invalid Speciationpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if (activeParts[i] == YES || nApplied == 0)
						{
						if (!strcmp(modelParams[i].speciationPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].speciationUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].speciationUni[0] >= modelParams[i].speciationUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Speciationpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].speciationUni[0], modelParams[i].speciationUni[1]);
								else
									MrBayesPrint ("%s   Setting Speciationpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].speciationUni[0], modelParams[i].speciationUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].speciationPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].speciationExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Speciationpr to Exponential(%1.2lf)\n", spacer, modelParams[i].speciationExp);
							else
								MrBayesPrint ("%s   Setting Speciationpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].speciationExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].speciationPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (AreDoublesEqual(tempD, 0.0, ETA)==YES)
								{
								MrBayesPrint ("%s   Speciation rate cannot be fixed to 0.0\n", spacer);
								return (ERROR);
								}
							modelParams[i].speciationFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Speciationpr to Fixed(%1.2lf)\n", spacer, modelParams[i].speciationFix);
							else
								MrBayesPrint ("%s   Setting Speciationpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].speciationFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Extinctionpr (extinctionPr) ****************************************************/
		else if (!strcmp(parmName, "Extinctionpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].extinctionPr, tempStr);
					}
				else
					{
					MrBayesPrint ("%s   Invalid Extinctionpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if (activeParts[i] == YES || nApplied == 0)
						{
						if (!strcmp(modelParams[i].extinctionPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].extinctionUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].extinctionUni[0] >= modelParams[i].extinctionUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Extinctionpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].extinctionUni[0], modelParams[i].extinctionUni[1]);
								else
									MrBayesPrint ("%s   Setting Extinctionpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].extinctionUni[0], modelParams[i].extinctionUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].extinctionPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].extinctionExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Extinctionpr to Exponential(%1.2lf)\n", spacer, modelParams[i].extinctionExp);
							else
								MrBayesPrint ("%s   Setting Extinctionpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].extinctionExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].extinctionPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].extinctionFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Extinctionpr to Fixed(%1.2lf)\n", spacer, modelParams[i].extinctionFix);
							else
								MrBayesPrint ("%s   Setting Extinctionpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].extinctionFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Sampleprob (sampleProb) *****************************************************/
		else if (!strcmp(parmName, "Sampleprob"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				if (tempD <= 1.0 && tempD > 0.0)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0))
							{
							modelParams[i].sampleProb = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Sampleprob to %1.2lf\n", spacer, modelParams[i].sampleProb);
							else
								MrBayesPrint ("%s   Setting Sampleprob to %1.2lf for partition %d\n", spacer, modelParams[i].sampleProb, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Sampleprob argument (should be between 0 and 1)\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Treeheightpr (treeHeightPr) ****************************************************/
		else if (!strcmp(parmName, "Treeheightpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].treeHeightPr, tempStr);
					}
				else
					{
					MrBayesPrint ("%s   Invalid Treeheightpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if (activeParts[i] == YES || nApplied == 0)
						{
						if (!strcmp(modelParams[i].treeHeightPr,"Gamma"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].treeHeightGamma[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Treeheightpr to Gamma(%1.2lf,%1.2lf)\n", spacer, modelParams[i].treeHeightGamma[0], modelParams[i].treeHeightGamma[1]);
								else
									MrBayesPrint ("%s   Setting Treeheightpr to Gamma(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].treeHeightGamma[0], modelParams[i].treeHeightGamma[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].treeHeightPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].treeHeightExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Treeheightpr to Exponential(%1.2lf)\n", spacer, modelParams[i].treeHeightExp);
							else
								MrBayesPrint ("%s   Setting Treeheightpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].treeHeightExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/*prior for joint gene tree*/
		else if (!strcmp(parmName, "BEST"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
                                
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);				
			       jointGenePr = tempInt;
				if(jointGenePr == 1)
					MrBayesPrint ("%s   Setting BEST to YES \n", spacer);
				else
					MrBayesPrint ("%s   Setting BEST to NO \n", spacer);

	              	expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}                 
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}
		else if (!strcmp(parmName, "PoissonMean"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
                                
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);				
			        poissonMean = tempD;
				
				MrBayesPrint ("%s   Setting poissonMean to %f \n", spacer,poissonMean);

	              		expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}                 
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}
		else if (!strcmp(parmName, "propTemp"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
                                
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);				
			        propTemp = tempD;
				
				MrBayesPrint ("%s   Setting proportion of generations when the temperature reaches 1.0 to %f \n", spacer,propTemp);

	              		expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}                 
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}

		else if (!strcmp(parmName, "PopmuPr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
                                
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);				
			        spMupr = tempInt;
				if(spMupr == 1)
					MrBayesPrint ("%s   Setting non-clock species tree to YES \n", spacer);
				else
					MrBayesPrint ("%s   Setting non-clock species tree to NO \n", spacer);

	              	expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}                 
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}
		else if (!strcmp(parmName, "PopmuPrAlpha"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
                                
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);				
			        spMupralpha = tempInt;
				if(spMupr == 1)
					MrBayesPrint ("%s   The alpha of the prior distribution of mutations is \n", spacer,spMupralpha);
				else
					{
					MrBayesPrint ("%s   Setting non-clock species tree to YES \n", spacer);
					exit(-1);
					}
	              	expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}                 
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}

		 else if (!strcmp(parmName, "SpeciesTreePr"))
                        {
                        if (expecting == Expecting(EQUALSIGN))
                                expecting = Expecting(NUMBER);

                        else if (expecting == Expecting(NUMBER))
                                {
                                sscanf (tkn, "%d", &tempInt);
                                speciestreePr = tempInt;
                                if(speciestreePr == 1)
                                        MrBayesPrint ("%s   Setting SpeciesTreePr to Birth and Death \n", spacer);
                                else
                                        MrBayesPrint ("%s   Setting SpeciesTreePr to uniform \n", spacer);

                        expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                                }
                        else if (expecting == Expecting(COMMA))
                                {
                                expecting  = Expecting(NUMBER);
                                }
                        else
                                return (ERROR);
                        }
		/* prior for joint gene tree */               
		/* set GeneMupr (GeneMuPr) **************************************************************/
		else if (!strcmp(parmName, "GeneMuPr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].GeneMuPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid GeneMuPr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN))
						{
						if (!strcmp(modelParams[i].GeneMuPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].GeneMuUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].GeneMuUni[0] >= modelParams[i].GeneMuUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting GeneMuPr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].GeneMuUni[0], modelParams[i].GeneMuUni[1]);
								else
									MrBayesPrint ("%s   Setting GeneMuPr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].GeneMuUni[0], modelParams[i].GeneMuUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].GeneMuPr,"Gamma"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].GeneMuGamma[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
                                                                                                                else
							               {
                                                                                                                                if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting GeneMuPr to Gamma(%1.2lf,%1.2lf)\n", spacer, modelParams[i].GeneMuGamma[0], modelParams[i].GeneMuGamma[1]);
							                else
								MrBayesPrint ("%s   Setting GeneMuPr to Gamma(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].GeneMuGamma[0], modelParams[i].GeneMuGamma[1], i+1);
							                expecting  = Expecting(RIGHTPAR);
                                                                                                                                }
							}
						else if (!strcmp(modelParams[i].GeneMuPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (fabs(tempD)<ETA)
								{
								MrBayesPrint ("%s   GeneMu cannot be fixed to 0.0\n", spacer);
								return (ERROR);
								}
							modelParams[i].GeneMuFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting GeneMuPr to Fixed(%1.2lf)\n", spacer, modelParams[i].GeneMuFix);
							else
								MrBayesPrint ("%s   Setting GeneMuPr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].GeneMuFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);			
		}
		/* set Thetapr (thetaPr) **************************************************************/
		else if (!strcmp(parmName, "Thetapr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].thetaPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Thetapr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA || modelParams[i].dataType == PROTEIN) )
						{
						if (!strcmp(modelParams[i].thetaPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].thetaUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].thetaUni[0] >= modelParams[i].thetaUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Thetapr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].thetaUni[0], modelParams[i].thetaUni[1]);
								else
									MrBayesPrint ("%s   Setting Thetapr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].thetaUni[0], modelParams[i].thetaUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						/*prior for joint gene tree*/
						else if (!strcmp(modelParams[i].thetaPr,"Gamma") || !strcmp(modelParams[i].thetaPr,"Invgamma"))
                                                        {
                                                        sscanf (tkn, "%lf", &tempD);
                                                        modelParams[i].thetaGamma[numVars[i]++] = tempD;
                                                        if (numVars[i] == 1)
                                                                expecting  = Expecting(COMMA);
                                                        else
                                                                {
                                                                if (nApplied == 0 && numCurrentDivisions == 1)
									{
									if(!strcmp(modelParams[i].thetaPr,"Invgamma"))
                                                                        	MrBayesPrint ("%s   Setting Thetapr to InvGamma(%1.2lf,%1.2lf) (for BEST)\n", spacer, modelParams[i].thetaGamma[0], modelParams[i].thetaGamma[1]);
									else
										MrBayesPrint ("%s   Setting Thetapr to Gamma(%1.2lf,%1.2lf) \n", spacer, modelParams[i].thetaGamma[0], modelParams[i].thetaGamma[1]);
                                                               		}
								 else
									{
									if (!strcmp(modelParams[i].thetaPr,"Gamma"))
                                                                        	MrBayesPrint ("%s   Setting Thetapr to Gamma(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].thetaGamma[0], modelParams[i].thetaGamma[1], i+1);
									else
										if(i == 0)
										MrBayesPrint ("%s   Setting Thetapr to InvGamma(%1.2lf,%1.2lf) (for BEST)\n", spacer, modelParams[i].thetaGamma[0], modelParams[i].thetaGamma[1], i+1);
									}
                                                                expecting  = Expecting(RIGHTPAR);
                                                                }
                                                        }
						else if (!strcmp(modelParams[i].thetaPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].thetaExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Thetapr to Exponential(%1.2lf)\n", spacer, modelParams[i].thetaExp);
							else
								MrBayesPrint ("%s   Setting Thetapr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].thetaExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].thetaPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (AreDoublesEqual(tempD, 0.0, ETA)==YES)
								{
								MrBayesPrint ("%s   Theta cannot be fixed to 0.0\n", spacer);
								return (ERROR);
								}
							modelParams[i].thetaFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Thetapr to Fixed(%1.2lf)\n", spacer, modelParams[i].thetaFix);
							else
								MrBayesPrint ("%s   Setting Thetapr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].thetaFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		
		/* set clock waiting time prior (calWaitPr) *********************************************************/
		else if (!strcmp(parmName, "Calwaitpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].calWaitPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Calwaitpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if (activeParts[i] == YES || nApplied == 0)
						{
						if (!strcmp(modelParams[i].calWaitPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].calWaitUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].calWaitUni[0] >= modelParams[i].calWaitUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Calwaitpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].thetaUni[0], modelParams[i].thetaUni[1]);
								else
									MrBayesPrint ("%s   Setting Calwaitpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].thetaUni[0], modelParams[i].thetaUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].calWaitPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].calWaitExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Calwaitpr to Exponential(%1.2lf)\n", spacer, modelParams[i].calWaitExp);
							else
								MrBayesPrint ("%s   Setting Calwaitpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].calWaitExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].calWaitPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (AreDoublesEqual(tempD, 0.0,ETA)==YES)
								{
								MrBayesPrint ("%s   Clock waiting time cannot be fixed to 0.0\n", spacer);
								return (ERROR);
								}
							modelParams[i].calWaitFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Calwaitpr to Fixed(%1.2lf)\n", spacer, modelParams[i].thetaFix);
							else
								MrBayesPrint ("%s   Setting Calwaitpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].thetaFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Growthpr (growthPr) **************************************************************/
		else if (!strcmp(parmName, "Growthpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				isNegative = NO;
				}
			else if (expecting == Expecting(ALPHA))
				{
				isNegative = NO;
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							strcpy(modelParams[i].growthPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Growthpr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER) | Expecting(DASH);
				}
			else if (expecting == Expecting(DASH))
				{
				expecting  = Expecting(NUMBER);
				isNegative = YES;
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].growthPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (isNegative == YES)
								tempD *= -1.0;
							modelParams[i].growthUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].growthUni[0] >= modelParams[i].growthUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Growthpr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].growthUni[0], modelParams[i].growthUni[1]);
								else
									MrBayesPrint ("%s   Setting Growthpr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].growthUni[0], modelParams[i].growthUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].growthPr,"Normal"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (isNegative == YES)
								tempD *= -1.0;
							modelParams[i].growthNorm[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].growthNorm[1] < 0.0)
									{
									MrBayesPrint ("%s   Variance for normal distribution should be greater than zero\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Growthpr to Normal(%1.2lf,%1.2lf)\n", spacer, modelParams[i].growthNorm[0], modelParams[i].growthNorm[1]);
								else
									MrBayesPrint ("%s   Setting Growthpr to Normal(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].growthNorm[0], modelParams[i].growthNorm[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].growthPr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].growthExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Growthpr to Exponential(%1.2lf)\n", spacer, modelParams[i].growthExp);
							else
								MrBayesPrint ("%s   Setting Growthpr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].growthExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].thetaPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							if (isNegative == YES)
								tempD *= -1.0;
							modelParams[i].growthFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Growthpr to Fixed(%1.2lf)\n", spacer, modelParams[i].growthFix);
							else
								MrBayesPrint ("%s   Setting Growthpr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].growthFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				isNegative = NO;
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Aamodelpr (aaModelPr) **************************************************************/
		else if (!strcmp(parmName, "Aamodelpr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				foundAaSetting = foundExp = modelIsFixed = foundDash = NO;
				fromI = 0;
				for (i=0; i<10; i++)
					tempAaModelPrs[i] = 0.0;
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (foundAaSetting == NO)
					{
					if (IsArgValid(tkn, tempStr) == NO_ERROR)
						{
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == PROTEIN)
								{
								strcpy(modelParams[i].aaModelPr, tempStr);
								if (!strcmp(modelParams[i].aaModelPr, "Mixed"))
									{
									if (nApplied == 0 && numCurrentDivisions == 1)
										MrBayesPrint ("%s   Setting Aamodelpr to %s\n", spacer, modelParams[i].aaModelPr);
									else
										MrBayesPrint ("%s   Setting Aamodelpr to %s for partition %d\n", spacer, modelParams[i].aaModelPr, i+1);
									}
								}
							}
						}
					else
						{
						MrBayesPrint ("%s   Invalid Aamodelpr argument\n", spacer);
						return (ERROR);
						}
					foundAaSetting = YES;
					if (!strcmp(tempStr, "Fixed"))
						{
						modelIsFixed = YES;
						expecting = Expecting(LEFTPAR);
						}
					else
						{
						expecting = Expecting(LEFTPAR) | Expecting(PARAMETER) | Expecting(SEMICOLON);
						}
					}
				else
					{
					if (modelIsFixed == YES)
						{
						if (IsSame ("Poisson", tkn) == SAME      || IsSame ("Poisson", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Poisson");
						else if (IsSame ("Equalin", tkn) == SAME || IsSame ("Equalin", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Equalin");
						else if (IsSame ("Jones", tkn) == SAME   || IsSame ("Jones", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Jones");
						else if (IsSame ("Dayhoff", tkn) == SAME || IsSame ("Dayhoff", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Dayhoff");
						else if (IsSame ("Mtrev", tkn) == SAME   || IsSame ("Mtrev", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Mtrev");
						else if (IsSame ("Mtmam", tkn) == SAME   || IsSame ("Mtmam", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Mtmam");
						else if (IsSame ("Wag", tkn) == SAME     || IsSame ("Wag", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Wag");
						else if (IsSame ("Rtrev", tkn) == SAME   || IsSame ("Rtrev", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Rtrev");
						else if (IsSame ("Cprev", tkn) == SAME   || IsSame ("Cprev", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Cprev");
						else if (IsSame ("Vt", tkn) == SAME      || IsSame ("Vt", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Vt");
						else if (IsSame ("Blosum", tkn) == SAME  || IsSame ("Blosum", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Blosum");
						else if (IsSame ("Blossum", tkn) == SAME  || IsSame ("Blossum", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Blosum");
						else if (IsSame ("Gtr", tkn) == SAME     || IsSame ("Gtr", tkn) == CONSISTENT_WITH)
							strcpy (tempStr, "Gtr");
						else
							{
							MrBayesPrint ("%s   Invalid amino acid model\n", spacer);
							return (ERROR);
							}
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == PROTEIN)
								{
								if (!strcmp(modelParams[i].aaModelPr, "Fixed"))
									{
									strcpy(modelParams[i].aaModel, tempStr);
									if (nApplied == 0 && numCurrentDivisions == 1)
										MrBayesPrint ("%s   Setting Aamodelpr to Fixed(%s)\n", spacer, modelParams[i].aaModel);
									else
										MrBayesPrint ("%s   Setting Aamodelpr to Fixed(%s) for partition %d\n", spacer, modelParams[i].aaModel, i+1);
									}
								else
									{
									MrBayesPrint ("%s   You cannot assign an amino acid matrix for mixed models\n", spacer);
									return (ERROR);
									}
								}
							}
						expecting = Expecting(RIGHTPAR);
						}
					else
						{
						if (IsSame ("Exponential", tkn) == SAME || IsSame ("Exponential", tkn) == CONSISTENT_WITH)
							{
							foundExp = YES;
							expecting = Expecting(LEFTPAR);
							}
						else	
							{
							MrBayesPrint ("%s   Invalid argument \"%s\"\n", spacer, tkn);
							return (ERROR);
							}
						}

					}
				}
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				if (fromI >= 10)
					{
					MrBayesPrint ("%s   Too many arguments in Aamodelpr\n", spacer);
					return (ERROR);
					}
				if (modelIsFixed == NO)
					{
					if (foundExp == YES)
						{
						if (foundDash == YES)
							tempAaModelPrs[fromI++] = -tempD;
						else
							tempAaModelPrs[fromI++] = tempD;
						expecting  = Expecting(RIGHTPAR);
						}
					else
						{
						if (foundDash == YES)
							{
							MrBayesPrint ("%s   Unexpected \"-\" in Aamodelpr\n", spacer);
							return (ERROR);
							}
						else
							{
							if (tempD <= 0.000000000001)
								tempAaModelPrs[fromI++] = -1000000000;
							else
								tempAaModelPrs[fromI++] = (MrBFlt) log(tempD);
							}
						expecting  = Expecting(COMMA) | Expecting(RIGHTPAR);
						}
					foundDash = NO;
					}
				else
					{
					MrBayesPrint ("%s   Not expecting a number\n", spacer);
					return (ERROR);
					}
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				if (modelIsFixed == YES)
					expecting  = Expecting(ALPHA);
				else
					{
					if (foundExp == YES)
						expecting  = Expecting(NUMBER) | Expecting(DASH);
					else
						expecting  = Expecting(NUMBER) | Expecting(ALPHA);
					}
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				if (modelIsFixed == YES)
					expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				else
					{
					if (foundExp == YES)
						{
						foundExp = NO;
						expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
						}
					else	
						{
						if (fromI < 10)
							{
							MrBayesPrint ("%s   Too few arguments in Aamodelpr\n", spacer);
							return (ERROR);
							}
						nApplied = NumActiveParts ();
						for (i=0; i<numCurrentDivisions; i++)
							{
							if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == PROTEIN)
								{
								if (!strcmp(modelParams[i].aaModelPr, "Fixed"))
									{
									MrBayesPrint ("%s   You cannot assign model prior probabilities for a fixed amino acid model\n", spacer);
									return (ERROR);
									}
								else
									{
									for (j=0; j<10; j++)
										modelParams[i].aaModelPrProbs[j] = tempAaModelPrs[j];
									}
								}
							}
						expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
						}
					}				
				}
			else if (expecting == Expecting(DASH))
				{
				if (foundExp == YES)
					{
					foundDash = YES;
					expecting  = Expecting(NUMBER);
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument \"%s\"\n", spacer, tkn);
					return (ERROR);
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				if (modelIsFixed == YES)
					{
					MrBayesPrint ("%s   Not expecting \"%s\"\n", spacer, tkn);
					return (ERROR);
					}
				else
					expecting  = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else
				return (ERROR);
			}		
		/* set Brownscalepr (brownScalesPr) ****************************************************/
		else if (!strcmp(parmName, "Brownscalepr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
							strcpy(modelParams[i].brownScalesPr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid Brownscalepr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && modelParams[i].dataType == CONTINUOUS)
						{
						if (!strcmp(modelParams[i].brownScalesPr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].brownScalesUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].brownScalesUni[0] >= modelParams[i].brownScalesUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Brownscalepr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].brownScalesUni[0], modelParams[i].brownScalesUni[1]);
								else
									MrBayesPrint ("%s   Setting Brownscalepr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].brownScalesUni[0], modelParams[i].brownScalesUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].brownScalesPr,"Gamma"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].brownScalesGamma[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting Brownscalepr to Gamma Mean=%1.2lf Var=%1.2lf\n", spacer, modelParams[i].brownScalesGamma[0], modelParams[i].brownScalesGamma[1]);
								else
									MrBayesPrint ("%s   Setting Brownscalepr to Gamma Mean=%1.2lf Var=%1.2lf for partition %d\n", spacer, modelParams[i].brownScalesGamma[0], modelParams[i].brownScalesGamma[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].brownScalesPr,"Gammamean"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].brownScalesGammaMean = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Brownscalepr to Gamma Mean=<char. ave.> Var=%1.2lf\n", spacer, modelParams[i].brownScalesGammaMean);
							else
								MrBayesPrint ("%s   Setting Brownscalepr to Gamma Mean=<char.ave.> Var=%1.2lf for partition %d\n", spacer, modelParams[i].brownScalesGammaMean, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].brownScalesPr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].brownScalesFix = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting Brownscalepr to Fixed(%1.2lf)\n", spacer, modelParams[i].brownScalesFix);
							else
								MrBayesPrint ("%s   Setting Brownscalepr to Fixed(%1.2lf) for partition %d\n", spacer, modelParams[i].brownScalesFix, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set M10betapr (m10betapr) ********************************************************/
		else if (!strcmp(parmName, "M10betapr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].m10betapr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid M10betapr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].m10betapr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m10betaUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].m10betaUni[0] >= modelParams[i].m10betaUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting M10betapr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].m10betaUni[0], modelParams[i].m10betaUni[1]);
								else
									MrBayesPrint ("%s   Setting M10betapr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].m10betaUni[0], modelParams[i].m10betaUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].m10betapr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m10betaExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting M10betapr to Exponential(%1.2lf)\n", spacer, modelParams[i].m10betaExp);
							else
								MrBayesPrint ("%s   Setting M10betapr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].m10betaExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].m10betapr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m10betaFix[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting M10betapr to Fixed(%1.2lf,%1.2lf)\n", spacer, modelParams[i].m10betaFix[0], modelParams[i].m10betaFix[1]);
								else
									MrBayesPrint ("%s   Setting M10betapr to Fixed(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].m10betaFix[0], modelParams[i].m10betaFix[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set M10gammapr (m10gammapr) ********************************************************/
		else if (!strcmp(parmName, "M10gammapr"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							strcpy(modelParams[i].m10gammapr, tempStr);
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid M10gammapr argument\n", spacer);
					return (ERROR);
					}
				expecting  = Expecting(LEFTPAR);
				for (i=0; i<MAX_NUM_DIVS; i++)
					numVars[i] = 0;
				}
			else if (expecting == Expecting(LEFTPAR))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(NUMBER))
				{
				nApplied = NumActiveParts ();
				for (i=0; i<numCurrentDivisions; i++)
					{
					if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
						{
						if (!strcmp(modelParams[i].m10gammapr,"Uniform"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m10gammaUni[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (modelParams[i].m10gammaUni[0] >= modelParams[i].m10gammaUni[1])
									{
									MrBayesPrint ("%s   Lower value for uniform should be greater than upper value\n", spacer);
									return (ERROR);
									}
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting M10gammapr to Uniform(%1.2lf,%1.2lf)\n", spacer, modelParams[i].m10gammaUni[0], modelParams[i].m10gammaUni[1]);
								else
									MrBayesPrint ("%s   Setting M10gammapr to Uniform(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].m10gammaUni[0], modelParams[i].m10gammaUni[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						else if (!strcmp(modelParams[i].m10gammapr,"Exponential"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m10gammaExp = tempD;
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting M10gammapr to Exponential(%1.2lf)\n", spacer, modelParams[i].m10gammaExp);
							else
								MrBayesPrint ("%s   Setting M10gammapr to Exponential(%1.2lf) for partition %d\n", spacer, modelParams[i].m10gammaExp, i+1);
							expecting  = Expecting(RIGHTPAR);
							}
						else if (!strcmp(modelParams[i].m10gammapr,"Fixed"))
							{
							sscanf (tkn, "%lf", &tempD);
							modelParams[i].m10gammaFix[numVars[i]++] = tempD;
							if (numVars[i] == 1)
								expecting  = Expecting(COMMA);
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Setting M10gammapr to Fixed(%1.2lf,%1.2lf)\n", spacer, modelParams[i].m10gammaFix[0], modelParams[i].m10gammaFix[1]);
								else
									MrBayesPrint ("%s   Setting M10gammapr to Fixed(%1.2lf,%1.2lf) for partition %d\n", spacer, modelParams[i].m10gammaFix[0], modelParams[i].m10gammaFix[1], i+1);
								expecting  = Expecting(RIGHTPAR);
								}
							}
						}
					}
				}
			else if (expecting == Expecting(COMMA))
				{
				expecting  = Expecting(NUMBER);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}



		
		
		
			

		else
			return (ERROR);
		}

	return (NO_ERROR);

}





int DoReport (void)

{

	return (NO_ERROR);
	
}





int DoReportParm (char *parmName, char *tkn)

{
	

	int			i, tempInt, nApplied;
	char		tempStr[100];

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before the report settings can be altered\n", spacer);
		return (ERROR);
		}
	if (inValidCommand == YES)
		{
		for (i=0; i<MAX_NUM_DIVS; i++)
			activeParts[i] = NO;
		inValidCommand = NO;
		}

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		/* set Applyto (Applyto) *************************************************************/
		if (!strcmp(parmName, "Applyto"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(LEFTPAR);
			else if (expecting == Expecting(LEFTPAR))
				{
				for (i=0; i<MAX_NUM_DIVS; i++)
					activeParts[i] = NO;
				fromI = toJ = -1;
				foundDash = NO;
				expecting = Expecting(NUMBER) | Expecting(ALPHA);
				}
			else if (expecting == Expecting(RIGHTPAR))
				{
				if (fromI != -1)
					activeParts[fromI-1] = YES;
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else if (expecting == Expecting(COMMA))
				{
				foundComma = YES;
				expecting = Expecting(NUMBER);
				}
			else if (expecting == Expecting(ALPHA))
				{
				if (IsSame ("All", tkn) == DIFFERENT)
					{
					MrBayesPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
					return (ERROR);
					}
				for (i=0; i<numCurrentDivisions; i++)
					activeParts[i] = YES;
				expecting  = Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempInt);
				if (tempInt > numCurrentDivisions)
					{
					MrBayesPrint ("%s   Partition delimiter is too large\n", spacer);
					return (ERROR);
					}
				if (fromI == -1)
					fromI = tempInt;
				else if (fromI != -1 && toJ == -1 && foundDash == YES && foundComma == NO)
					{
					toJ = tempInt;
					for (i=fromI-1; i<toJ; i++)
						activeParts[i] = YES;
					fromI = toJ = -1;
					foundDash = NO;
					}
				else if (fromI != -1 && toJ == -1 && foundDash == NO && foundComma == YES)
					{
					activeParts[fromI-1] = YES;
					fromI = tempInt;
					foundComma = NO;
					}
					
				expecting  = Expecting(COMMA);
				expecting |= Expecting(DASH);
				expecting |= Expecting(RIGHTPAR);
				}
			else if (expecting == Expecting(DASH))
				{
				foundDash = YES;
				expecting = Expecting(NUMBER);
				}
			else
				return (ERROR);
			}
		/* set report format of tratio ***************************************************/
		else if (!strcmp(parmName, "Tratio"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting (ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					tempInt = NO;
					for (i=0; i<numCurrentDivisions; i++)
						{
						/* check that data type is correct; we do not know yet if the user will specify
						a nst=2 model so we cannot check that tratio is an active parameter */
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].tratioFormat, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting transition/transversion rate ratio (tratio) format to %s\n", spacer, modelParams[i].tratioFormat);
							else
								MrBayesPrint ("%s   Setting transition/transversion rate ratio (tratio) format to %s for partition %d\n", spacer, modelParams[i].tratioFormat, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid transition/transversion rate ratio (tratio) format \n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set report format of revmat ***************************************************/
		else if (!strcmp(parmName, "Revmat"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting (ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					tempInt = NO;
					for (i=0; i<numCurrentDivisions; i++)
						{
						/* check that data type is correct; we do not know yet if the user will specify
						   a nst=6 model so we cannot check that revmat is an active parameter */
						if ((activeParts[i] == YES || nApplied == 0) && (modelParams[i].dataType == DNA || modelParams[i].dataType == RNA))
							{
							strcpy(modelParams[i].revmatFormat, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting reversible rate matrix (revmat) format to %s\n", spacer, modelParams[i].revmatFormat);
							else
								MrBayesPrint ("%s   Setting reversible rate matrix (revmat) format to %s for partition %d\n", spacer, modelParams[i].revmatFormat, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid reversible rate matrix (revmat) format \n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set report format of ratemult *************************************************/
		else if (!strcmp(parmName, "Ratemult"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting (ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					tempInt = NO;
					for (i=0; i<numCurrentDivisions; i++)
						{
						/* we do not know yet if the user will specify variable rates across partitions 
						   so only check that we have more than one partition in the model */
						if ((activeParts[i] == YES || nApplied == 0) && numCurrentDivisions > 1)
							{
							strcpy(modelParams[i].ratemultFormat, tempStr);
							if (nApplied == 0 && numCurrentDivisions == 1)
								MrBayesPrint ("%s   Setting rate multiplier (ratemult) format to %s\n", spacer, modelParams[i].ratemultFormat);
							else
								MrBayesPrint ("%s   Setting rate multiplier (ratemult) format to %s for partition %d\n", spacer, modelParams[i].ratemultFormat, i+1);
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid rate multiplier (ratemult) format \n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set inferancstates ***********************************************************/
		else if (!strcmp(parmName, "Ancstates"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting (ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							{
							strcpy(modelParams[i].inferAncStates,tempStr);
							if (!strcmp(tempStr,"Yes"))
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Reporting ancestral states (if applicable)\n", spacer);
								else
									MrBayesPrint ("%s   Reporting ancestral states for partition %d (if applicable)\n", spacer, i+1);
								}
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Not reporting ancestral states\n", spacer);
								else
									MrBayesPrint ("%s   Not reporting ancestral states for partition %d\n", spacer, i+1);
								}
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid ancstates option\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set inferSiteRates ***************************************************************/
		else if (!strcmp(parmName, "Siterates"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting (ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							{
							strcpy (modelParams[i].inferSiteRates, tempStr);
							if (!strcmp(tempStr,"Yes"))
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Reporting site rates (if applicable)\n", spacer);
								else
									MrBayesPrint ("%s   Reporting site rates for partition %d (if applicable)\n", spacer, i+1);
								}
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Not reporting site rates\n", spacer);
								else
									MrBayesPrint ("%s   Not reporting site rates for partition %d\n", spacer, i+1);
								}
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid siterates option\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set inferpossel *************************************************************/
		else if (!strcmp(parmName, "Possel"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting (ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					nApplied = NumActiveParts ();
					for (i=0; i<numCurrentDivisions; i++)
						{
						if (activeParts[i] == YES || nApplied == 0)
							{
							strcpy (modelParams[i].inferPosSel, tempStr);
							if (!strcmp(tempStr, "Yes"))
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Reporting positive selection (if applicable)\n", spacer);
								else
									MrBayesPrint ("%s   Reporting positive selection for partition %d (if applicable)\n", spacer, i+1);
								}
							else
								{
								if (nApplied == 0 && numCurrentDivisions == 1)
									MrBayesPrint ("%s   Not reporting positive selection\n", spacer);
								else
									MrBayesPrint ("%s   Not reporting positive selection for partition %d\n", spacer, i+1);
								}
							}
						}
					}
				else
					{
					MrBayesPrint ("%s   Invalid inferpossel option\n", spacer);
					return (ERROR);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		else
			{
			return (ERROR);
			}
		}

	return (NO_ERROR);		
}




int DoUnlink (void)

{

	int			i, j;
	
	MrBayesPrint ("%s   Unlinking\n", spacer);
	
	/* update status of linkTable */
	for (j=0; j<NUM_LINKED; j++)
		{
		for (i=0; i<numCurrentDivisions; i++)
			{
			if (tempLinkUnlink[j][i] == YES)
				{
				linkTable[j][i] = ++linkNum;
				}
			}
		}
	
#	if 0
	for (j=0; j<NUM_LINKED; j++)
		{
		MrBayesPrint ("%s   ", spacer);
		for (i=0; i<numCurrentDivisions; i++)
			MrBayesPrint ("%d", linkTable[j][i]);
		MrBayesPrint ("\n");
		}
#	endif

	/* reinitialize the temporary table */
	for (j=0; j<NUM_LINKED; j++)
		for (i=0; i<MAX_NUM_DIVS; i++)
			tempLinkUnlink[j][i] = NO;

	return (NO_ERROR);
	
}





int DoShowModel (void)

{

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before the model can be defined\n", spacer);
		return (ERROR);
		}

	if (CheckModel() == ERROR)
		return (ERROR);

	return (NO_ERROR);

}





int IsModelSame (int whichParam, int part1, int part2, int *isApplic1, int *isApplic2)

{

	int			i, isSame, isFirstNucleotide, isSecondNucleotide, isFirstProtein, isSecondProtein, nDiff, temp1, temp2;

	isSame = YES;
	*isApplic1 = YES;
	*isApplic2 = YES;
	isFirstNucleotide = isSecondNucleotide = NO;
	if (modelParams[part1].dataType == DNA || modelParams[part1].dataType == RNA)
		isFirstNucleotide = YES;
	if (modelParams[part2].dataType == DNA || modelParams[part2].dataType == RNA)
		isSecondNucleotide = YES;		
	isFirstProtein = isSecondProtein = NO;
	if (modelParams[part1].dataType == PROTEIN)
		isFirstProtein = YES;
	if (modelParams[part2].dataType == PROTEIN)
		isSecondProtein = YES;		
	
	if (whichParam == P_TRATIO)
		{
		/* Check the ti/tv rate ratio for partitions 1 and 2. */

		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and tratio does not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and tratio does not apply */
		
		/* Check that the data are nucleotide for both partitions 1 and 2 */
		if (isFirstNucleotide == NO)
			*isApplic1 = NO; /* part1 is not nucleotide data so tratio does not apply */
		if (isSecondNucleotide == NO)
			*isApplic2 = NO; /* part2 is not nucleotide data so tratio does not apply */
		
		/* check that nst=2 for both partitions */
		if (strcmp(modelParams[part1].nst, "2"))
			*isApplic1 = NO; /* part1 does not have nst=2 and tratio does not apply */
		if (strcmp(modelParams[part2].nst, "2"))
			*isApplic2 = NO; /* part2 does not have nst=2 and tratio does not apply */
		
		/* Check if part1 & part2 are restriction */
		if (modelParams[part1].dataType == RESTRICTION)
			*isApplic1 = NO;
		if (modelParams[part2].dataType == RESTRICTION)
			*isApplic2 = NO;

		/* If Nst = 2 for both part1 and part2, we now need to check if the prior is the same for both. */
		if (!strcmp(modelParams[part1].tRatioPr,"Beta") && !strcmp(modelParams[part2].tRatioPr,"Beta"))
			{
			if (AreDoublesEqual (modelParams[part1].tRatioDir[0], modelParams[part2].tRatioDir[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].tRatioDir[1], modelParams[part2].tRatioDir[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].tRatioPr,"Fixed") && !strcmp(modelParams[part2].tRatioPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].tRatioFix, modelParams[part2].tRatioFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check to see if tratio is inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if tratio is inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_REVMAT)
		{
		/* Check the GTR rates for partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and GTR rates do not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and GTR rates do not apply */

		/* Check that the data are nucleotide or protein for both partitions 1 and 2 */
		if (isFirstNucleotide == NO && isFirstProtein == NO)
			*isApplic1 = NO; /* part1 is not nucleotide or protein data so GTR rates do not apply */
		if (isSecondNucleotide == NO && isSecondProtein == NO)
			*isApplic2 = NO; /* part2 is not nucleotide or protein data so GTR rates do not apply */

		/* check that nst=6 for both partitions if nucleotide */
		if (isFirstNucleotide == YES && strcmp(modelParams[part1].nst, "6") && strcmp(modelParams[part1].nst, "Mixed"))
			*isApplic1 = NO; /* part1 does not have nst=6/Mixed and GTR rates do not apply */
		if (isSecondNucleotide == YES && strcmp(modelParams[part2].nst, "6") && strcmp(modelParams[part2].nst, "Mixed"))
			*isApplic2 = NO; /* part2 does not have nst=6/Mixed and GTR rates do not apply */
			
		/* check that model is GTR for both partitions if protein */
		if (isFirstProtein == YES && (strcmp(modelParams[part1].aaModel,"Gtr")!=0 || strcmp(modelParams[part1].aaModelPr,"Fixed")!=0))
			*isApplic1 = NO;
		if (isSecondProtein == YES && (strcmp(modelParams[part2].aaModel,"Gtr")!=0 || strcmp(modelParams[part2].aaModelPr,"Fixed")!=0))
			*isApplic2 = NO;

		/* check that data type is the same for both partitions */
		if (isFirstNucleotide == YES && isSecondNucleotide == NO)
			isSame = NO;
		if (isFirstProtein == YES && isSecondProtein == NO)
			isSame = NO;

		/* GTR applies to both part1 and part2. We now need to check if the prior is the same for both. */
		if (isFirstNucleotide == YES)
			{
			if (!strcmp(modelParams[part1].revMatPr,"Dirichlet") && !strcmp(modelParams[part2].revMatPr,"Dirichlet"))
				{
				for (i=0; i<6; i++)
					{
					if (AreDoublesEqual (modelParams[part1].revMatDir[i], modelParams[part2].revMatDir[i], (MrBFlt) 0.00001) == NO)
						isSame = NO;
					}
				}
			else if (!strcmp(modelParams[part1].revMatPr,"Fixed") && !strcmp(modelParams[part2].revMatPr,"Fixed"))
				{
				for (i=0; i<6; i++)
					{
					if (AreDoublesEqual (modelParams[part1].revMatFix[i], modelParams[part2].revMatFix[i], (MrBFlt) 0.00001) == NO)
						isSame = NO;
					}
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}
		else /* if (isFirstProtein == YES) */
			{
			if (!strcmp(modelParams[part1].aaRevMatPr,"Dirichlet") && !strcmp(modelParams[part2].aaRevMatPr,"Dirichlet"))
				{
				for (i=0; i<190; i++)
					{
					if (AreDoublesEqual (modelParams[part1].aaRevMatDir[i], modelParams[part2].aaRevMatDir[i], (MrBFlt) 0.00001) == NO)
						isSame = NO;
					}
				}
			else if (!strcmp(modelParams[part1].aaRevMatPr,"Fixed") && !strcmp(modelParams[part2].aaRevMatPr,"Fixed"))
				{
				for (i=0; i<190; i++)
					{
					if (AreDoublesEqual (modelParams[part1].aaRevMatFix[i], modelParams[part2].aaRevMatFix[i], (MrBFlt) 0.00001) == NO)
						isSame = NO;
					}
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}

		/* Check to see if the GTR rates are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if GTR rates are inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_OMEGA)
		{
		/* Check the nonsynonymous/synonymous rate ratio for partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and nonsynonymous/synonymous rate ratio does not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and nonsynonymous/synonymous rate ratio does not apply */
		
		/* Check that the data are nucleotide for both partitions 1 and 2 */
		if (isFirstNucleotide == NO)
			*isApplic1 = NO; /* part1 is not nucleotide data so a nonsynonymous/synonymous rate ratio does not apply */
		if (isSecondNucleotide == NO)
			*isApplic2 = NO; /* part2 is not nucleotide data so a nonsynonymous/synonymous rate ratio does not apply */
		
		/* Check that the model structure is the same for both. The nucmodel should be "codon". */
		if (strcmp(modelParams[part1].nucModel, "Codon"))
			*isApplic1 = NO; /* part1 does not have Nucmodel = Codon and nonsynonymous/synonymous rate ratio does not apply */
		if (strcmp(modelParams[part2].nucModel, "Codon"))
			*isApplic2 = NO; /* part2 does not have Nucmodel = Codon and nonsynonymous/synonymous rate ratio does not apply */
		
		/* Assuming that Nucmodel = Codon for both part1 and part2, we now need to check if the prior is the
		   same for both. */
		if (!strcmp(modelParams[part1].omegaVar, "M3") && !strcmp(modelParams[part2].omegaVar, "M3"))
			{
			if (!strcmp(modelParams[part1].m3omegapr, "Exponential") && !strcmp(modelParams[part2].m3omegapr, "Exponential"))
				{
				}
			else if (!strcmp(modelParams[part1].m3omegapr, "Fixed") && !strcmp(modelParams[part1].m3omegapr, "Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].m3omegaFixed[0], modelParams[part2].m3omegaFixed[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].m3omegaFixed[1], modelParams[part2].m3omegaFixed[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

			if (!strcmp(modelParams[part1].codonCatFreqPr, "Dirichlet") && !strcmp(modelParams[part2].codonCatFreqPr, "Dirichlet"))
				{
				if (AreDoublesEqual (modelParams[part1].codonCatDir[0], modelParams[part2].codonCatDir[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatDir[1], modelParams[part2].codonCatDir[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatDir[2], modelParams[part2].codonCatDir[2], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].codonCatFreqPr, "Fixed") && !strcmp(modelParams[part1].codonCatFreqPr, "Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[0], modelParams[part2].codonCatFreqFix[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[1], modelParams[part2].codonCatFreqFix[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[2], modelParams[part2].codonCatFreqFix[2], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}
		else if (!strcmp(modelParams[part1].omegaVar, "M10") && !strcmp(modelParams[part2].omegaVar, "M10"))
			{			
			if (!strcmp(modelParams[part1].codonCatFreqPr, "Dirichlet") && !strcmp(modelParams[part2].codonCatFreqPr, "Dirichlet"))
				{
				if (AreDoublesEqual (modelParams[part1].codonCatDir[0], modelParams[part2].codonCatDir[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatDir[1], modelParams[part2].codonCatDir[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].codonCatFreqPr, "Fixed") && !strcmp(modelParams[part1].codonCatFreqPr, "Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[0], modelParams[part2].codonCatFreqFix[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[1], modelParams[part2].codonCatFreqFix[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}
		else if (!strcmp(modelParams[part1].omegaVar, "Ny98") && !strcmp(modelParams[part2].omegaVar, "Ny98"))
			{
			if (!strcmp(modelParams[part1].ny98omega1pr, "Beta") && !strcmp(modelParams[part2].ny98omega1pr, "Beta"))
				{
				if (AreDoublesEqual (modelParams[part1].ny98omega1Beta[0], modelParams[part2].ny98omega1Beta[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].ny98omega1Beta[1], modelParams[part2].ny98omega1Beta[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].ny98omega1pr, "Fixed") && !strcmp(modelParams[part1].ny98omega1pr, "Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].ny98omega1Fixed, modelParams[part2].ny98omega1Fixed, (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			if (!strcmp(modelParams[part1].ny98omega3pr, "Uniform") && !strcmp(modelParams[part2].ny98omega3pr, "Uniform"))
				{
				if (AreDoublesEqual (modelParams[part1].ny98omega3Uni[0], modelParams[part2].ny98omega3Uni[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].ny98omega3Uni[1], modelParams[part2].ny98omega3Uni[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].ny98omega3pr, "Exponential") && !strcmp(modelParams[part1].ny98omega3pr, "Exponential"))
				{
				if (AreDoublesEqual (modelParams[part1].ny98omega3Exp, modelParams[part2].ny98omega3Exp, (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].ny98omega3pr, "Fixed") && !strcmp(modelParams[part1].ny98omega3pr, "Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].ny98omega3Fixed, modelParams[part2].ny98omega3Fixed, (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			if (!strcmp(modelParams[part1].codonCatFreqPr, "Dirichlet") && !strcmp(modelParams[part2].codonCatFreqPr, "Dirichlet"))
				{
				if (AreDoublesEqual (modelParams[part1].codonCatDir[0], modelParams[part2].codonCatDir[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatDir[1], modelParams[part2].codonCatDir[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatDir[2], modelParams[part2].codonCatDir[2], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].codonCatFreqPr, "Fixed") && !strcmp(modelParams[part1].codonCatFreqPr, "Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[0], modelParams[part2].codonCatFreqFix[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[1], modelParams[part2].codonCatFreqFix[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].codonCatFreqFix[2], modelParams[part2].codonCatFreqFix[2], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}
		else if (!strcmp(modelParams[part1].omegaVar, "Equal") && !strcmp(modelParams[part2].omegaVar, "Equal"))
			{
			if (!strcmp(modelParams[part1].omegaPr,"Dirichlet") && !strcmp(modelParams[part2].omegaPr,"Dirichlet"))
				{
				if (AreDoublesEqual (modelParams[part1].omegaDir[0], modelParams[part2].omegaDir[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].omegaDir[1], modelParams[part2].omegaDir[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else if (!strcmp(modelParams[part1].omegaPr,"Fixed") && !strcmp(modelParams[part2].omegaPr,"Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].omegaFix, modelParams[part2].omegaFix, (MrBFlt) 0.00001) == NO)
					isSame = NO;
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
		
		/* Check to see if omega is inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if omega is inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_PI)
		{
		/* Check the state frequencies for partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and state frequencies do not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and state frequencies do not apply */

		/* Check that the data are not CONTINUOUS for partitions 1 and 2 */
		if (modelParams[part1].dataType == CONTINUOUS)
			*isApplic1 = NO; /* state frequencies do not make sense for part1 */
		if (modelParams[part2].dataType == CONTINUOUS)
			*isApplic2 = NO; /* state frequencies do not make sense for part2 */
			
		/* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
		if (isFirstNucleotide != isSecondNucleotide)
			isSame = NO; /* data are not both nucleotide or both note nucleotide */
		else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
			isSame = NO; /* data are not the same */

		/* Check that the model structure is the same for both partitions */
		if (strcmp(modelParams[part1].nucModel, modelParams[part2].nucModel))
			isSame = NO; /* the nucleotide models are different */
		if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel) && !(!strcmp(modelParams[part1].nucModel, "Codon") && !strcmp(modelParams[part2].nucModel, "Codon")))
			isSame = NO; /* the models have different covarion struture */
			
		/* If both partitions have nucmodel=codon, then we also have to make certain that the same genetic code is used. */
		if (!strcmp(modelParams[part1].nucModel, "Codon") && !strcmp(modelParams[part2].nucModel, "Codon"))
			{
			if (strcmp(modelParams[part1].geneticCode,modelParams[part2].geneticCode))
				isSame = NO; /* the models have different genetic codes */
			}
		
		/* Let's see if the prior is the same. */
		if (modelParams[part1].dataType == STANDARD && modelParams[part2].dataType == STANDARD)
			{
			/* The data are morphological (STANDARD). The state frequencies are specified by a
			   symmetric beta distribution, the parameter of which needs to be the same to apply to both
			   partitions. Note that symPiPr = -1 is equivalent to setting the variance to 0.0. */
			if (!strcmp(modelParams[part1].symPiPr,"Uniform") && !strcmp(modelParams[part2].symPiPr,"Uniform"))
				{
				if (AreDoublesEqual (modelParams[part1].symBetaUni[0], modelParams[part2].symBetaUni[0], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].symBetaUni[1], modelParams[part2].symBetaUni[1], (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (modelParams[part1].numBetaCats != modelParams[part2].numBetaCats)
					isSame = NO;	/* can't link because the discrete beta approximation is different */
				}
			else if (!strcmp(modelParams[part1].symPiPr,"Exponential") && !strcmp(modelParams[part2].symPiPr,"Exponential"))
				{
				if (AreDoublesEqual (modelParams[part1].symBetaExp, modelParams[part2].symBetaExp, (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (modelParams[part1].numBetaCats != modelParams[part2].numBetaCats)
					isSame = NO;	/* can't link because the discrete beta approximation is different */
				}
			else if (!strcmp(modelParams[part1].symPiPr,"Fixed") && !strcmp(modelParams[part2].symPiPr,"Fixed"))
				{
				if (AreDoublesEqual (modelParams[part1].symBetaFix, modelParams[part2].symBetaFix, (MrBFlt) 0.00001) == NO)
					isSame = NO;
				if (AreDoublesEqual (modelParams[part1].symBetaFix, (MrBFlt) -1.0, (MrBFlt) 0.00001) == NO && modelParams[part1].numBetaCats != modelParams[part2].numBetaCats)
					isSame = NO;	/* can't link because the discrete beta approximation is different */
				}
			else
				isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
			}
		if (modelParams[part1].dataType == PROTEIN && modelParams[part2].dataType == PROTEIN)
			{
			/* We are dealing with protein data. */
			if (!strcmp(modelParams[part1].aaModelPr, modelParams[part2].aaModelPr))
				{
				if (!strcmp(modelParams[part1].aaModelPr, "Fixed"))
					{
					/* only have a single, fixed, amino acid rate matrix */
					if (!strcmp(modelParams[part1].aaModel, modelParams[part2].aaModel))
						{}
					else
						isSame = NO; /* we have different amino acid models, and the state frequencies must be different */
					/* if we have an equalin model or Gtr model, then we need to check the prior on the state frequencies */
					if (!strcmp(modelParams[part1].aaModel, "Equalin") || !strcmp(modelParams[part1].aaModel, "Gtr"))
						{
						if (!strcmp(modelParams[part1].stateFreqPr, modelParams[part2].stateFreqPr))
							{
							/* the prior form is the same */
							if (!strcmp(modelParams[part1].stateFreqPr, "Dirichlet")) /* both prior models must be dirichlet */
								{
								for (i=0; i<modelParams[part1].nStates; i++)
									if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (MrBFlt) 0.00001) == NO)
										isSame = NO; /* the dirichlet parameters are different */
								}
							else /* both prior models must be fixed */
								{
								if (!strcmp(modelParams[part1].stateFreqsFixType, modelParams[part2].stateFreqsFixType))
									{
									/* if (!strcmp(modelParams[part1].stateFreqsFixType, "Empirical"))
										isSame = NO;     Even though it is unlikely that the empirical values for both partitions are exactly the same, we will
										                 allow isSame to equal YES. This means pooled base frequencies are used to determine the empirical
										                 base frequencies. The user can still unlink this parameter. */
									if (!strcmp(modelParams[part1].stateFreqsFixType, "User"))
										{
										for (i=0; i<modelParams[part1].nStates; i++)
											if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (MrBFlt) 0.00001) == NO)
												isSame = NO; /* the user-specified base frequencies are different */
										}
									/* if the frequencies were both fixed to "equal", we are golden, and they are the same */
									}
								else
									isSame = NO; /* the fixed parameters must be the same. The only other possibility is that the
									                user specified equal or empirical for one partition and then specified specific
									                numbers (user) for the other _and_ happened to set the user values to the equal
									                or empirical values. We ignore this possibility. */
								}
							}
						}
					}
				else
					{
					/* averaging over models */
					if (linkTable[P_AAMODEL][part1] != linkTable[P_AAMODEL][part2])
						isSame = NO; /* the amino acid model is mixed, but independently estimated */
					}
				}
			}
		else
			{
			/* Otherwise, we are dealing with RESTRICTION or NUCLEOTIDE data. The dirichlet should be the same
			   for both partitions. */
			if (!strcmp(modelParams[part1].stateFreqPr, modelParams[part2].stateFreqPr))
				{
				/* the prior form is the same */
				if (!strcmp(modelParams[part1].stateFreqPr, "Dirichlet")) /* both prior models must be dirichlet */
					{
					for (i=0; i<modelParams[part1].nStates; i++)
						if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (MrBFlt) 0.00001) == NO)
							isSame = NO; /* the dirichlet parameters are different */
					}
				else /* both prior models must be fixed */
					{
					if (!strcmp(modelParams[part1].stateFreqsFixType, modelParams[part2].stateFreqsFixType))
						{
						/* if (!strcmp(modelParams[part1].stateFreqsFixType, "Empirical"))
							isSame = NO;     Even though it is unlikely that the empirical values for both partitions are exactly the same, we will
							                 allow isSame to equal YES. This means pooled base frequencies are used to determine the empirical
							                 base frequencies. The user can still unlink this parameter. */
						if (!strcmp(modelParams[part1].stateFreqsFixType, "User"))
							{
							for (i=0; i<modelParams[part1].nStates; i++)
								if (AreDoublesEqual (modelParams[part1].stateFreqsDir[i], modelParams[part2].stateFreqsDir[i], (MrBFlt) 0.00001) == NO)
									isSame = NO; /* the user-specified base frequencies are different */
							}
						/* if the frequencies were both fixed to "equal", we are golden, and they are the same */
						}
					else
						isSame = NO; /* the fixed parameters must be the same. The only other possibility is that the
						                user specified equal or empirical for one partition and then specified specific
						                numbers (user) for the other _and_ happened to set the user values to the equal
						                or empirical values. We ignore this possibility. */
					}
				}
			}

		/* Check to see if the state frequencies are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the state frequencies are inapplicable for either partition, then the parameter cannot be the same */
		
		}
	else if (whichParam == P_SHAPE)
		{
		/* Check the gamma shape parameter for partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and gamma shape parameter does not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and gamma shape parameter does not apply */

		/* Check that the data are not CONTINUOUS for partitions 1 and 2 */
		if (modelParams[part1].dataType == CONTINUOUS)
			*isApplic1 = NO; /* the gamma shape parameter does not make sense for part1 */
		if (modelParams[part2].dataType == CONTINUOUS)
			*isApplic2 = NO; /* the gamma shape parameter does not make sense for part2 */

		/* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
		if (isFirstNucleotide != isSecondNucleotide)
			isSame = NO; /* data are not both nucleotide */
		else if ( modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
			isSame = NO; /* data are not the same */

		/* Let's check that the gamma shape parameter is even relevant for the two partitions */
		if (!strcmp(modelParams[part1].ratesModel, "Equal") || !strcmp(modelParams[part1].ratesModel, "Propinv"))
			*isApplic1 = NO; /* the gamma shape parameter does not make sense for part1 */
		if (!strcmp(modelParams[part2].ratesModel, "Equal") || !strcmp(modelParams[part2].ratesModel, "Propinv"))
			*isApplic2 = NO; /* the gamma shape parameter does not make sense for part2 */
		
		/* We may have a nucleotide model. Make certain the models are not of type codon. */
		if (!strcmp(modelParams[part1].nucModel, "Codon"))
			*isApplic1 = NO; /* we have a codon model for part1, and a gamma shape parameter does not apply */
		if (!strcmp(modelParams[part2].nucModel, "Codon"))
			*isApplic2 = NO;/* we have a codon model for part2, and a gamma shape parameter does not apply */

		/* Check that the model structure is the same for both partitions */
		if ((!strcmp(modelParams[part1].nucModel, "4by4") || !strcmp(modelParams[part1].nucModel, "Doublet")) && !strcmp(modelParams[part2].nucModel, "Codon"))
			isSame = NO; /* the nucleotide models are incompatible with the same shape parameter */
		if ((!strcmp(modelParams[part2].nucModel, "4by4") || !strcmp(modelParams[part2].nucModel, "Doublet")) && !strcmp(modelParams[part1].nucModel, "Codon"))
			isSame = NO; /* the nucleotide models are incompatible with the same shape parameter */
		/*if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
			isSame = NO;*/ /* the models have different covarion struture */  /* NOTE: Perhaps we should allow the possiblity that
			                                                                         the gamma shape parameter is the same for the
			                                                                         case where one partition has a covarion model
			                                                                         but the other does not and both datatypes are
			                                                                         the same. */
		
		/* Check that the number of rate categories is the same */
		if (modelParams[part1].numGammaCats != modelParams[part2].numGammaCats)
			isSame = NO; /* the number of rate categories is not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check that the priors are the same. */
		if (!strcmp(modelParams[part1].shapePr,"Uniform") && !strcmp(modelParams[part2].shapePr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].shapeUni[0], modelParams[part2].shapeUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].shapeUni[1], modelParams[part2].shapeUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].shapePr,"Exponential") && !strcmp(modelParams[part2].shapePr,"Exponential"))
			{
			if (AreDoublesEqual (modelParams[part1].shapeExp, modelParams[part2].shapeExp, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].shapePr,"Fixed") && !strcmp(modelParams[part2].shapePr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].shapeFix, modelParams[part2].shapeFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
		
		/* Check to see if the shape parameter is inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the shape parameter is inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_PINVAR)
		{
		/* Check the proportion of invariable sites parameter for partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and proportion of invariable sites parameter does not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and proportion of invariable sites parameter does not apply */

		/* Check that the data are not CONTINUOUS for partitions 1 and 2 */
		if (modelParams[part1].dataType == CONTINUOUS)
			*isApplic1 = NO; /* the proportion of invariable sites parameter does not make sense for part1 */
		if (modelParams[part2].dataType == CONTINUOUS)
			*isApplic2 = NO; /* the proportion of invariable sites parameter does not make sense for part2 */

		/* Now, check that the data are the same (i.e., both nucleotide or both amino acid, or whatever). */
		if (isFirstNucleotide != isSecondNucleotide)
			isSame = NO; /* data are not both nucleotide */
		else if ( modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
			isSame = NO; /* data are not the same */

		/* Let's check that proportion of invariable sites parameter is even relevant for the two partitions */
		if (!strcmp(modelParams[part1].ratesModel, "Equal") || !strcmp(modelParams[part1].ratesModel, "Gamma") || !strcmp(modelParams[part1].ratesModel, "Adgamma"))
			*isApplic1 = NO; /* the proportion of invariable sites parameter does not make sense for part1 */
		if (!strcmp(modelParams[part2].ratesModel, "Equal") || !strcmp(modelParams[part2].ratesModel, "Gamma") || !strcmp(modelParams[part2].ratesModel, "Adgamma"))
			*isApplic2 = NO; /* the proportion of invariable sites parameter does not make sense for part2 */
			
		/* It is not sensible to have a covarion model and a proportion of invariable sites */
		if (!strcmp(modelParams[part1].covarionModel, "Yes"))
			*isApplic1 = NO;
		if (!strcmp(modelParams[part2].covarionModel, "Yes"))
			*isApplic2 = NO;
		
		/* We have a nucleotide model. Make certain the models are not of type codon. */
		if (!strcmp(modelParams[part1].nucModel, "Codon"))
			*isApplic1 = NO; /* we have a codon model for part1, and a proportion of invariable sites parameter does not apply */
		if (!strcmp(modelParams[part2].nucModel, "Codon"))
			*isApplic2 = NO;/* we have a codon model for part2, and a proportion of invariable sites parameter does not apply */

		/* Check that the model structure is the same for both partitions */
		if (strcmp(modelParams[part1].nucModel, modelParams[part2].nucModel))
			isSame = NO; /* the nucleotide models are different */
		if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
			isSame = NO; /* the models have different covarion struture */
		
		/* check the priors */
		if (!strcmp(modelParams[part1].pInvarPr,"Uniform") && !strcmp(modelParams[part2].pInvarPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].pInvarUni[0], modelParams[part2].pInvarUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].pInvarUni[1], modelParams[part2].pInvarUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].pInvarPr,"Fixed") && !strcmp(modelParams[part2].pInvarPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].pInvarFix, modelParams[part2].pInvarFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */
		
		/* Check to see if the switching rates are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the switching rates are inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_CORREL)
		{
		/* Check the autocorrelation parameter for gamma rates on partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and autocorrelation parameter does not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and autocorrelation parameter does not apply */

		/* Check that the data are either DNA, RNA, or PROTEIN for partitions 1 and 2 */
		if (modelParams[part1].dataType != DNA && modelParams[part1].dataType != RNA && modelParams[part1].dataType != PROTEIN)
			*isApplic1 = NO; /* the switching rates do not make sense for part1 */
		if (modelParams[part2].dataType != DNA && modelParams[part2].dataType != RNA && modelParams[part2].dataType != PROTEIN)
			*isApplic2 = NO; /* the switching rates do not make sense for part2 */
			
		/* Now, check that the data are the same (i.e., both nucleotide or both amino acid). */
		if (isFirstNucleotide != isSecondNucleotide)
			isSame = NO; /* one or the other is nucleotide, so they cannot be the same */
		else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
			isSame = NO; /* data are not both nucleotide or both amino acid */

		/* Let's check that autocorrelation parameter is even relevant for the two partitions */
		if (strcmp(modelParams[part1].ratesModel, "Adgamma"))
			*isApplic1 = NO; /* the autocorrelation parameter does not make sense for part1 */
		if (strcmp(modelParams[part2].ratesModel, "Adgamma"))
			*isApplic2 = NO; /* the autocorrelation parameter does not make sense for part2 */

		/* Assuming that we have a nucleotide model, make certain the models are not of type codon. */
		if (!strcmp(modelParams[part1].nucModel, "Codon"))
			*isApplic1 = NO; /* we have a codon model for part1, and a autocorrelation parameter does not apply */
		if (!strcmp(modelParams[part2].nucModel, "Codon"))
			*isApplic2 = NO; /* we have a codon model for part2, and a autocorrelation parameter does not apply */
		
		/* Check that the model structure is the same for both partitions */
		if (strcmp(modelParams[part1].nucModel, modelParams[part2].nucModel))
			isSame = NO; /* the nucleotide models are different */
		if (strcmp(modelParams[part1].covarionModel, modelParams[part2].covarionModel))
			isSame = NO; /* the models have different covarion struture */

		/* Check the priors for both partitions. */
		if (!strcmp(modelParams[part1].adGammaCorPr,"Uniform") && !strcmp(modelParams[part2].adGammaCorPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].corrUni[0], modelParams[part2].corrUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].corrUni[1], modelParams[part2].corrUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].adGammaCorPr,"Fixed") && !strcmp(modelParams[part2].adGammaCorPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].corrFix, modelParams[part2].corrFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check to see if the switching rates are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the switching rates are inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_SWITCH)
		{
		/* Check the covarion switching rates on partitions 1 and 2. */
		
		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and switching rates do not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and switching rates do not apply */
		
		/* Check that the data are either DNA, RNA, or PROTEIN for partitions 1 and 2 */
		if (modelParams[part1].dataType != DNA && modelParams[part1].dataType != RNA && modelParams[part1].dataType != PROTEIN)
			*isApplic1 = NO; /* the switching rates do not make sense for part1 */
		if (modelParams[part2].dataType != DNA && modelParams[part2].dataType != RNA && modelParams[part2].dataType != PROTEIN)
			*isApplic2 = NO; /* the switching rates do not make sense for part2 */
			
		/* Now, check that the data are the same (i.e., both nucleotide or both amino acid). */
		if (isFirstNucleotide != isSecondNucleotide)
			isSame = NO; /* one or the other is nucleotide, so they cannot be the same */
		else if (modelParams[part1].dataType != modelParams[part2].dataType && isFirstNucleotide == NO)
			isSame = NO; /* data are not both nucleotide or both amino acid */

		/* Lets check that covarion model has been selected for partitions 1 and 2 */
		if (!strcmp(modelParams[part1].covarionModel, "No"))
			*isApplic1 = NO; /* the switching rates do not make sense for part1 */
		if (!strcmp(modelParams[part2].covarionModel, "No"))
			*isApplic2 = NO; /* the switching rates do not make sense for part2 */

		/* If we have a nucleotide model make certain the models are not of type codon or doublet. */
		if (!strcmp(modelParams[part1].nucModel, "Codon") || !strcmp(modelParams[part1].nucModel, "Doublet"))
			*isApplic1 = NO; /* we have a codon model for part1, and a autocorrelation parameter does not apply */
		if (!strcmp(modelParams[part2].nucModel, "Codon") || !strcmp(modelParams[part2].nucModel, "Doublet"))
			*isApplic2 = NO; /* we have a codon model for part2, and a autocorrelation parameter does not apply */

		/* Check that the priors are the same. */
		if (!strcmp(modelParams[part1].covSwitchPr,"Uniform") && !strcmp(modelParams[part2].covSwitchPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].covswitchUni[0], modelParams[part2].covswitchUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].covswitchUni[1], modelParams[part2].covswitchUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].covSwitchPr,"Exponential") && !strcmp(modelParams[part2].covSwitchPr,"Exponential"))
			{
			if (AreDoublesEqual (modelParams[part1].covswitchExp, modelParams[part2].covswitchExp, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].covSwitchPr,"Fixed") && !strcmp(modelParams[part2].covSwitchPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].covswitchFix[0], modelParams[part2].covswitchFix[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].covswitchFix[1], modelParams[part2].covswitchFix[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check to see if the switching rates are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the switching rates are inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_RATEMULT)
		{
		/* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
		   silly and doesn't take this information into account) and a rate multiplier is nonsensical. */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; 
		
		/* Check that the branch lengths are at least proportional. */
		if (IsModelSame (P_BRLENS, part1, part2, &temp1, &temp2) == NO)
			isSame = NO;
		if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
			isSame = NO;

		/* See if the rate prior is the same for the partitions */
		if (strcmp(modelParams[part1].ratePr, modelParams[part2].ratePr) != 0)
			isSame = NO;

		/* See if either rate is fixed to 1.0 */
		if (!strcmp(modelParams[part1].ratePr, "Fixed"))
			*isApplic1 = NO;
		if (!strcmp(modelParams[part2].ratePr, "Fixed"))
			*isApplic2 = NO;
					
		/* Now, check that there is more than one partition. In SetModel, we call this function with both part1 and part2
		   the same. If part1 = part2, then we know that we have only one division, and that a rate multiplier is not
		   relevant. */
		if (part1 == part2)
			*isApplic1 = *isApplic2 = NO;
			
		/* Check to see if rate multipliers are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; 
			
		}
	else if (whichParam == P_TOPOLOGY)
		{
		/* Check the topology for partitions 1 and 2. */
		
		/* If the prior is different, then the topologies cannot be the same. */
		if (strcmp(modelParams[part1].topologyPr, modelParams[part2].topologyPr))
			isSame = NO;

		/* If both partitions have topologies constrained, then we need to make certain that the constraints are the same. */
		if (!strcmp(modelParams[part1].topologyPr, "Constraints") && !strcmp(modelParams[part2].topologyPr, "Constraints"))
			{
			if (modelParams[part1].numActiveConstraints != modelParams[part2].numActiveConstraints)
				isSame = NO;
			else
				{
				nDiff = 0;
				for (i=0; i<modelParams[part1].numActiveConstraints; i++)
					if (modelParams[part1].activeConstraints[i] != modelParams[part2].activeConstraints[i])
						nDiff++;
				if (nDiff != 0)
					isSame = NO;
				}
			}

		/* Check if one model is parsimony but the other is not, then the topologies cannot be the same */
		/*if (strcmp(modelParams[part1].parsModel, modelParams[part2].parsModel))
			isSame = NO; */ /* I turned this off, because, on thinking about it, I see no reason why we cannot have the same topology
			when one partition is parsimony and the other is not. However, in this case, we cannot have the same branch lengths,
			so we worry about that later. */
			
		if (!strcmp(modelParams[part1].parsModel, "No") && !strcmp(modelParams[part2].parsModel, "No"))
			{
			/* We are dealing with a real substitution models for both partitions */
			
			/* The branch length prior should be the same */
			if (strcmp(modelParams[part1].brlensPr, modelParams[part2].brlensPr))
				isSame = NO;
				
			/* if both partitions have are unconstrained brlens, then we need to check that the priors on the branch lengths are the same */
			if (!strcmp(modelParams[part1].brlensPr, "Unconstrained") && !strcmp(modelParams[part2].brlensPr, "Unconstrained"))
				{
				if (strcmp(modelParams[part1].unconstrainedPr, modelParams[part2].unconstrainedPr))
					isSame = NO;
				else
					{
					if (!strcmp(modelParams[part1].unconstrainedPr, "Uniform"))
						{
						if (AreDoublesEqual (modelParams[part1].brlensUni[0], modelParams[part2].brlensUni[0], (MrBFlt) 0.00001) == NO)
							isSame = NO;
						if (AreDoublesEqual (modelParams[part1].brlensUni[1], modelParams[part2].brlensUni[1], (MrBFlt) 0.00001) == NO)
							isSame = NO;
						}
					else
						{
						if (AreDoublesEqual (modelParams[part1].brlensExp, modelParams[part2].brlensExp, (MrBFlt) 0.00001) == NO)
							isSame = NO;
						}
					}
				}
			
			/* if both partitions have clock brlens, then we need to check that the priors on the clock are the same */
			if (!strcmp(modelParams[part1].brlensPr, "Clock") && !strcmp(modelParams[part2].brlensPr, "Clock"))
				{
				if (strcmp(modelParams[part1].clockPr, modelParams[part2].clockPr))
					isSame = NO;
				else
					{
					if (!strcmp(modelParams[part1].clockPr, "Birthdeath"))
						{
						if (!strcmp(modelParams[part1].speciationPr,"Uniform") && !strcmp(modelParams[part2].speciationPr,"Uniform"))
							{
							if (AreDoublesEqual (modelParams[part1].speciationUni[0], modelParams[part2].speciationUni[0], (MrBFlt) 0.00001) == NO)
								isSame = NO;
							if (AreDoublesEqual (modelParams[part1].speciationUni[1], modelParams[part2].speciationUni[1], (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else if (!strcmp(modelParams[part1].speciationPr,"Exponential") && !strcmp(modelParams[part2].speciationPr,"Exponential"))
							{
							if (AreDoublesEqual (modelParams[part1].speciationExp, modelParams[part2].speciationExp, (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else if (!strcmp(modelParams[part1].speciationPr,"Fixed") && !strcmp(modelParams[part2].speciationPr,"Fixed"))
							{
							if (AreDoublesEqual (modelParams[part1].speciationFix, modelParams[part2].speciationFix, (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else
							isSame = NO;

						if (!strcmp(modelParams[part1].extinctionPr,"Uniform") && !strcmp(modelParams[part2].extinctionPr,"Uniform"))
							{
							if (AreDoublesEqual (modelParams[part1].extinctionUni[0], modelParams[part2].extinctionUni[0], (MrBFlt) 0.00001) == NO)
								isSame = NO;
							if (AreDoublesEqual (modelParams[part1].extinctionUni[1], modelParams[part2].extinctionUni[1], (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else if (!strcmp(modelParams[part1].extinctionPr,"Exponential") && !strcmp(modelParams[part2].extinctionPr,"Exponential"))
							{
							if (AreDoublesEqual (modelParams[part1].extinctionExp, modelParams[part2].extinctionExp, (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else if (!strcmp(modelParams[part1].extinctionPr,"Fixed") && !strcmp(modelParams[part2].extinctionPr,"Fixed"))
							{
							if (AreDoublesEqual (modelParams[part1].extinctionFix, modelParams[part2].extinctionFix, (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else
							isSame = NO;
						}
					else if (!strcmp(modelParams[part1].clockPr, "Coalescence"))
						{
						if (!strcmp(modelParams[part1].thetaPr,"Uniform") && !strcmp(modelParams[part2].thetaPr,"Uniform"))
							{
							if (AreDoublesEqual (modelParams[part1].thetaUni[0], modelParams[part2].thetaUni[0], (MrBFlt) 0.00001) == NO)
								isSame = NO;
							if (AreDoublesEqual (modelParams[part1].thetaUni[1], modelParams[part2].thetaUni[1], (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else if (!strcmp(modelParams[part1].thetaPr,"Exponential") && !strcmp(modelParams[part2].thetaPr,"Exponential"))
							{
							if (AreDoublesEqual (modelParams[part1].thetaExp, modelParams[part2].thetaExp, (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else if (!strcmp(modelParams[part1].thetaPr,"Fixed") && !strcmp(modelParams[part2].thetaPr,"Fixed"))
							{
							if (AreDoublesEqual (modelParams[part1].thetaFix, modelParams[part2].thetaFix, (MrBFlt) 0.00001) == NO)
								isSame = NO;
							}
						else
							isSame = NO;
						}
					}
				}



			}
		}
	/* prior for joint gene tree */
	else if (whichParam == P_GENEMU)                                
		{
			if ( jointGenePr == 0 && strcmp(modelParams[part1].GeneMuPr,"Fixed"))		                     
			{ 
				MrBayesPrint ("%s   GeneMuPr should be fixed if BEST is not 1\n", spacer);
				exit(-1);
			}
			/*else if(spMupr==1 && strcmp(modelParams[part1].GeneMuPr,"Fixed"))		                     
			{ 
				MrBayesPrint ("%s   GeneMuPr should be fixed when non-clock species tree popmupr=1\n", spacer);
				exit(-1);
			}*/
			else
			{ 
				if (!strcmp(modelParams[part1].GeneMuPr,"Uniform") && !strcmp(modelParams[part2].GeneMuPr,"Uniform"))
				{
					if (AreDoublesEqual (modelParams[part1].GeneMuUni[0], modelParams[part2].GeneMuUni[0], (MrBFlt) 0.00001) == NO)
						isSame = NO;
					if (AreDoublesEqual (modelParams[part1].GeneMuUni[1], modelParams[part2].GeneMuUni[1], (MrBFlt) 0.00001) == NO)
						isSame = NO;
				}
				else if (!strcmp(modelParams[part1].GeneMuPr,"Gamma") && !strcmp(modelParams[part2].GeneMuPr,"Gamma"))
				{
					if (AreDoublesEqual (modelParams[part1].GeneMuGamma[0], modelParams[part2].GeneMuGamma[0], (MrBFlt) 0.00001) == NO)
						isSame = NO;
					if (AreDoublesEqual (modelParams[part1].GeneMuGamma[1], modelParams[part2].GeneMuGamma[1], (MrBFlt) 0.00001) == NO)
						isSame = NO;
				}
				else if (!strcmp(modelParams[part1].GeneMuPr,"Fixed") && !strcmp(modelParams[part2].GeneMuPr,"Fixed"))
				{
					if (AreDoublesEqual (modelParams[part1].GeneMuFix, modelParams[part2].GeneMuFix, (MrBFlt) 0.00001) == NO)
						isSame = NO;
				}
				else
					isSame = NO;
			}
		}
                              
	else if (whichParam == P_BRLENS)
		{
		/* Check the branch lengths for partitions 1 and 2. */

		/* First, if the topologies are different, the same branch lengths cannot apply. Note that the topologies will be different
		   if they differ in their prior information, such as the range on the speciation and extinction rates. This single check
		   does a lot. */
		if (IsModelSame (P_TOPOLOGY, part1, part2, &temp1, &temp2) == NO)
			isSame = NO;
		if (linkTable[P_TOPOLOGY][part1] != linkTable[P_TOPOLOGY][part2])
			isSame = NO;

		/* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
		   silly and doesn't take this information into account). */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO;

		/* Check to see if the branch lengths are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; 

		}
	else if (whichParam == P_SPECRATE)
		{
		/* Check the speciation rates for partitions 1 and 2. */

		/* First, if the branch lengths are different, we don't want to apply the same speciation rate to it */
		/*if (IsModelSame (P_BRLENS, part1, part2, &temp1, &temp2) == NO)
			isSame = NO;
		if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
			isSame = NO;*/ /* I think that it may make some sense, after all, to allow different branch lengths in
		  different partitions, even if the topology/branch lengths are different. */

		/* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
		   silly and doesn't take this information into account) and a speciation rate cannot be estimated. */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; 
			
		/* Check that the branch length prior is a clock:birthdeath for both partitions. */
		if (strcmp(modelParams[part1].brlensPr, "Clock"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].brlensPr, "Clock"))
			*isApplic2 = NO;
		if (strcmp(modelParams[part1].clockPr, "Birthdeath"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].clockPr, "Birthdeath"))
			*isApplic2 = NO;
		
		/* Now, check that the prior on the speciation rates are the same. */
		if (!strcmp(modelParams[part1].speciationPr,"Uniform") && !strcmp(modelParams[part2].speciationPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].speciationUni[0], modelParams[part2].speciationUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].speciationUni[1], modelParams[part2].speciationUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].speciationPr,"Exponential") && !strcmp(modelParams[part2].speciationPr,"Exponential"))
			{
			if (AreDoublesEqual (modelParams[part1].speciationExp, modelParams[part2].speciationExp, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].speciationPr,"Fixed") && !strcmp(modelParams[part2].speciationPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].speciationFix, modelParams[part2].speciationFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO;
		
		/* Check to see if the speciation rates are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO;
		/*prior for joint gene tree*/
		if(jointGenePr == 1 && speciestreePr == 1)
			isSame = YES; 

		}
	else if (whichParam == P_EXTRATE)
		{
		/* Check the extinction rates for partitions 1 and 2. */

		/* First, if the branch lengths are different, we don't want to apply the same extinction rate to it */
		/*if (IsModelSame (P_BRLENS, part1, part2, &temp1, &temp2) == NO)
			isSame = NO;
		if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
			isSame = NO;*/ /* I think that it may make some sense, after all, to allow different branch lengths in
		  different partitions, even if the topology/branch lengths are different. */

		/* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
		   silly and doesn't take this information into account) and a extinction rate cannot be estimated. */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; 
			
		/* Check that the branch length prior is a clock:birthdeath for both partitions. */
		if (strcmp(modelParams[part1].brlensPr, "Clock"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].brlensPr, "Clock"))
			*isApplic2 = NO;
		if (strcmp(modelParams[part1].clockPr, "Birthdeath"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].clockPr, "Birthdeath"))
			*isApplic2 = NO;
		
		/* Now, check that the prior on the extinction rates are the same. */
		if (!strcmp(modelParams[part1].extinctionPr,"Uniform") && !strcmp(modelParams[part2].extinctionPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].extinctionUni[0], modelParams[part2].extinctionUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].extinctionUni[1], modelParams[part2].extinctionUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].extinctionPr,"Exponential") && !strcmp(modelParams[part2].extinctionPr,"Exponential"))
			{
			if (AreDoublesEqual (modelParams[part1].extinctionExp, modelParams[part2].extinctionExp, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].extinctionPr,"Fixed") && !strcmp(modelParams[part2].extinctionPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].extinctionFix, modelParams[part2].extinctionFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO;
		
		/* Check to see if the speciation rates are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; 

		/*prior for joint gene tree*/
                if(jointGenePr == 1 && speciestreePr == 1)
                        isSame = YES;

		}
	else if (whichParam == P_THETA)
		{
		/* Check theta for partitions 1 and 2. */

		/* First, if the branch lengths are different, we don't want to apply the same theta to it */
		/*if (IsModelSame (P_BRLENS, part1, part2, &temp1, &temp2) == NO)
			isSame = NO;
		if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
			isSame = NO;*/ /* I think that it may make some sense, after all, to allow different branch lengths in
		  different partitions, even if the topology/branch lengths are different. */

		/* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
		   silly and doesn't take this information into account) and theta cannot be estimated. */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; 
			
		/* Check that the branch length prior is a clock:coalescence for both partitions. */
		if (strcmp(modelParams[part1].brlensPr, "Clock"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].brlensPr, "Clock"))
			*isApplic2 = NO;
		if (strcmp(modelParams[part1].clockPr, "Coalescence"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].clockPr, "Coalescence"))
			*isApplic2 = NO;
		
		/* Now, check that the prior on theta is the same. */
		if (!strcmp(modelParams[part1].thetaPr,"Uniform") && !strcmp(modelParams[part2].thetaPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].thetaUni[0], modelParams[part2].thetaUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].thetaUni[1], modelParams[part2].thetaUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].thetaPr,"Exponential") && !strcmp(modelParams[part2].thetaPr,"Exponential"))
			{
			if (AreDoublesEqual (modelParams[part1].thetaExp, modelParams[part2].thetaExp, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].thetaPr,"Fixed") && !strcmp(modelParams[part2].thetaPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].thetaFix, modelParams[part2].thetaFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO;
		
		/* Check to see if theta is inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; 

		}
	else if (whichParam == P_GROWTH)
		{
		/* Check growth rate for partitions 1 and 2. */

		/* First, if the branch lengths are different, we don't want to apply the same growth rate to it */
		/*if (IsModelSame (P_BRLENS, part1, part2, &temp1, &temp2) == NO)
			isSame = NO;
		if (linkTable[P_BRLENS][part1] != linkTable[P_BRLENS][part2])
			isSame = NO;*/ /* I think that it may make some sense, after all, to allow different branch lengths in
		  different partitions, even if the topology/branch lengths are different. */

		/* Check if the model is parsimony for either partition. If so, then the branch lengths cannot apply (as parsimony is very
		   silly and doesn't take this information into account) and growth rate cannot be estimated. */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; 
			
		/* Check that the branch length prior is a clock:coalescence for both partitions. */
		if (strcmp(modelParams[part1].brlensPr, "Clock"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].brlensPr, "Clock"))
			*isApplic2 = NO;
		if (strcmp(modelParams[part1].clockPr, "Coalescence"))
			*isApplic1 = NO;
		if (strcmp(modelParams[part2].clockPr, "Coalescence"))
			*isApplic2 = NO;
		
		/* Now, check that the prior on growth rate is the same. */
		if (!strcmp(modelParams[part1].growthPr,"Uniform") && !strcmp(modelParams[part2].growthPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].growthUni[0], modelParams[part2].growthUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].growthUni[1], modelParams[part2].growthUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].growthPr,"Exponential") && !strcmp(modelParams[part2].growthPr,"Exponential"))
			{
			if (AreDoublesEqual (modelParams[part1].growthExp, modelParams[part2].growthExp, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].growthPr,"Fixed") && !strcmp(modelParams[part2].growthPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].growthFix, modelParams[part2].growthFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO;
		
		/* Check to see if growth rate is inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; 

		}
	else if (whichParam == P_AAMODEL)
		{
		/* Check the amino acid model settings for partitions 1 and 2. */

		/* Check if the model is parsimony for either partition */
		if (!strcmp(modelParams[part1].parsModel, "Yes"))
			*isApplic1 = NO; /* part1 has a parsimony model and aamodel does not apply */
		if (!strcmp(modelParams[part2].parsModel, "Yes"))
			*isApplic2 = NO; /* part2 has a parsimony model and aamodel does not apply */
		
		/* Check that the data are protein for both partitions 1 and 2 */
		if (isFirstProtein == NO)
			*isApplic1 = NO; /* part1 is not amino acid data so tratio does not apply */
		if (isSecondProtein == NO)
			*isApplic2 = NO; /* part2 is not amino acid data so tratio does not apply */
			
		/* If the model is fixed for a partition, then it is not a free parameter and
		   we set it to isApplic = NO */
		if (!strcmp(modelParams[part1].aaModelPr,"Fixed"))
			*isApplic1 = NO; 
		if (!strcmp(modelParams[part2].aaModelPr,"Fixed"))
			*isApplic2 = NO; 

		/* We now need to check if the prior is the same for both. */
		if (!strcmp(modelParams[part1].aaModelPr,"Mixed") && !strcmp(modelParams[part2].aaModelPr,"Mixed"))
			{
			}
		else if (!strcmp(modelParams[part1].aaModelPr,"Fixed") && !strcmp(modelParams[part2].aaModelPr,"Fixed"))
			{
			if (strcmp(modelParams[part1].aaModel,modelParams[part2].aaModel))
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check to see if amino acid model is inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if tratio is inapplicable for either partition, then the parameter cannot be the same */

		}
	else if (whichParam == P_BRCORR)
		{
		/* Check the correlation parameter for brownian motion 1 and 2. */
		
		/* Check that the data are either CONTINUOUS for partitions 1 and 2 */
		if (modelParams[part1].dataType != CONTINUOUS)
			*isApplic1 = NO; /* the correlation parameter does not make sense for part1 */
		if (modelParams[part2].dataType != CONTINUOUS)
			*isApplic2 = NO; /* the correlation parameter does not make sense for part2 */
			
		/* Now, check that the data are the same. */
		if (modelParams[part1].dataType != modelParams[part2].dataType)
			isSame = NO; /* data are not both continuous */

		/* Check the priors for both partitions. */
		if (!strcmp(modelParams[part1].brownCorPr,"Uniform") && !strcmp(modelParams[part2].brownCorPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].brownCorrUni[0], modelParams[part2].brownCorrUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].brownCorrUni[1], modelParams[part2].brownCorrUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].brownCorPr,"Fixed") && !strcmp(modelParams[part2].brownCorPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].brownCorrFix, modelParams[part2].brownCorrFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check to see if the correlation parameters are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the correlation parameters are inapplicable for either partition, then the parameter cannot be the same */
		}
	else if (whichParam == P_BRSIGMA)
		{
		/* Check the sigma parameter for brownian motion 1 and 2. */
		
		/* Check that the data are either CONTINUOUS for partitions 1 and 2 */
		if (modelParams[part1].dataType != CONTINUOUS)
			*isApplic1 = NO; /* the sigma parameter does not make sense for part1 */
		if (modelParams[part2].dataType != CONTINUOUS)
			*isApplic2 = NO; /* the sigma parameter does not make sense for part2 */
			
		/* Now, check that the data are the same. */
		if (modelParams[part1].dataType != modelParams[part2].dataType)
			isSame = NO; /* data are not both continuous */

		/* Check the priors for both partitions. */
		if (!strcmp(modelParams[part1].brownScalesPr,"Uniform") && !strcmp(modelParams[part2].brownScalesPr,"Uniform"))
			{
			if (AreDoublesEqual (modelParams[part1].brownScalesUni[0], modelParams[part2].brownScalesUni[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].brownScalesUni[1], modelParams[part2].brownScalesUni[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].brownScalesPr,"Fixed") && !strcmp(modelParams[part2].brownScalesPr,"Fixed"))
			{
			if (AreDoublesEqual (modelParams[part1].brownScalesFix, modelParams[part2].brownScalesFix, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].brownScalesPr,"Gamma") && !strcmp(modelParams[part2].brownScalesPr,"Gamma"))
			{
			if (AreDoublesEqual (modelParams[part1].brownScalesGamma[0], modelParams[part2].brownScalesGamma[0], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			if (AreDoublesEqual (modelParams[part1].brownScalesGamma[1], modelParams[part2].brownScalesGamma[1], (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else if (!strcmp(modelParams[part1].brownScalesPr,"Gammamean") && !strcmp(modelParams[part2].brownScalesPr,"Gammamean"))
			{
			if (AreDoublesEqual (modelParams[part1].brownScalesGammaMean, modelParams[part2].brownScalesGammaMean, (MrBFlt) 0.00001) == NO)
				isSame = NO;
			}
		else
			isSame = NO; /* the priors are not the same, so we cannot set the parameter to be equal for both partitions */

		/* Check to see if the sigma parameters are inapplicable for either partition. */
		if ((*isApplic1) == NO || (*isApplic2) == NO)
			isSame = NO; /* if the sigma parameters are inapplicable for either partition, then the parameter cannot be the same */
		}
	else
		{
		MrBayesPrint ("%s   Could not find parameter in IsModelSame\n", spacer);
		return (NO);
		}
	
	return (isSame);
	
}





int InitializeLinks (void)

{

	int			i, j;
	
	linkNum = 0;
	for (i=0; i<NUM_LINKED; i++)
		{
		for (j=0; j<MAX_NUM_DIVS; j++)
			linkTable[i][j] = linkNum;
		}

	return (NO_ERROR);
	
}





int Link (void)

{

	int			i, j;
	
	for (j=0; j<NUM_LINKED; j++)
		{
		MrBayesPrint ("%4d -- ", j+1);
		for (i=0; i<numCurrentDivisions; i++)
			MrBayesPrint (" %2d", tempLinkUnlink[j][i]);
		MrBayesPrint ("\n");
		}
		
	return (NO_ERROR);
	
}





int NumActiveParts (void)

{

	int		i, nApplied;
	
	nApplied = 0;
	for (i=0; i<numCurrentDivisions; i++)
		if (activeParts[i] == YES)
			nApplied++;

	return (nApplied);
	
}





int NumStates (int part)

{

	if (modelParams[part].dataType == DNA || modelParams[part].dataType == RNA)
		{
		if (!strcmp(modelParams[part].nucModel, "4by4"))
			return (4);
		else if (!strcmp(modelParams[part].nucModel, "Doublet"))
			return (16);
		else
			{
			if (!strcmp(modelParams[part].geneticCode, "Universal"))
				return (61);
			else if (!strcmp(modelParams[part].geneticCode, "Vertmt"))
				return (60);
			else if (!strcmp(modelParams[part].geneticCode, "Mycoplasma"))
				return (62);
			else if (!strcmp(modelParams[part].geneticCode, "Yeast"))
				return (62);
			else if (!strcmp(modelParams[part].geneticCode, "Ciliates"))
				return (63);
			else if (!strcmp(modelParams[part].geneticCode, "Metmt"))
				return (62);
			}
		}
	else if (modelParams[part].dataType == PROTEIN)
		{
		return (20);
		}
	else if (modelParams[part].dataType == RESTRICTION)
		{
		return (2);
		}
	else if (modelParams[part].dataType == STANDARD)
		{
		return (10);
		}
		
	return (-1);

}





int SetUpModel (void)

{

	int 		i, j, k, m, modelId[MAX_NUM_DIVS], paramCount, isApplicable1, isApplicable2, isFirst, isSame, check[MAX_NUM_DIVS];

	for (j=0; j<NUM_LINKED; j++)
		for (i=0; i<numCurrentDivisions; i++)
			activeParams[j][i] = 0;
	
	if (numCurrentDivisions > 1)
		{
		paramCount = 0;
		for (j=0; j<NUM_LINKED; j++) /* loop over parameters */
			{
			isFirst = YES;
			for (i=0; i<numCurrentDivisions; i++)
				modelId[i] = 0;		
			for (i=0; i<numCurrentDivisions-1; i++) /* loop over partitions */
				{
				for (k=i+1; k<numCurrentDivisions; k++)
					{
					if (IsModelSame (j, i, k, &isApplicable1, &isApplicable2) == NO || linkTable[j][i] != linkTable[j][k])
						{
						/* we cannot link the parameters */
						if (isApplicable1 == NO)
							modelId[i] = -1;
						if (isApplicable2 == NO)
							modelId[k] = -1;
						if (isApplicable1 == YES)
							{
							if (isFirst == YES && modelId[i] == 0)
								{
								modelId[i] = ++paramCount;
								isFirst = NO;
								}
							else
								{
								if (modelId[i] == 0)
									modelId[i] = ++paramCount;
								}
							}
						if (modelId[k] == 0 && isApplicable2 == YES)
							modelId[k] = ++paramCount;
						}
					else
						{
						/* we can link the parameters */
						if (isFirst == YES)
							{
							if (modelId[i] == 0)
								modelId[i] = ++paramCount;
							isFirst = NO;
							}
						else
							{
							if (modelId[i] == 0)
								modelId[i] = ++paramCount;
							}
						modelId[k] = modelId[i];
						}
					}
				}
			for (i=0; i<numCurrentDivisions; i++)
				activeParams[j][i] = modelId[i];
			}
		}
	else
		{
		/* if we have only one partition, then we do things a bit differently */
		paramCount = 0;
		for (j=0; j<NUM_LINKED; j++) /* loop over parameters */
			{
			IsModelSame (j, 0, 0, &isApplicable1, &isApplicable2);
			if (isApplicable1 == YES)
				activeParams[j][0] = ++paramCount;
			else
				activeParams[j][0] = -1;
			}
		}
		
	/* This might be a good place to deal with the rate multiplier and some other weird things the user might
	   specify regarding topologies. */
	for (i=0; i<numCurrentDivisions; i++)
		{
		m = activeParams[P_RATEMULT][i];
		k = 0;
		for (j=0; j<numCurrentDivisions; j++)
			if (activeParams[P_RATEMULT][j] == m)
				k++;
		if (k == 1)
			activeParams[P_RATEMULT][i] = -1;
		}
	   
	/* Check that the same report format is specified for all partitions with the same rate multiplier */
	for (i=0; i<numCurrentDivisions; i++)
		check[i] = NO;
	for (i=0; i<numCurrentDivisions; i++)
		{
		m = activeParams[P_RATEMULT][i];
		if (m == -1 || check[i] == YES)
			continue;
		isSame = YES;
		for (j=i+1; j<numCurrentDivisions; j++)
			{
			if (activeParams[P_RATEMULT][j] == m)
				{
				check[i] = YES;
				if (strcmp(modelParams[i].ratemultFormat,modelParams[j].ratemultFormat)!= 0)
					{
					isSame = NO;
					strcpy (modelParams[j].ratemultFormat, modelParams[i].ratemultFormat);
					}
				}
			}
		if (isSame == NO)
			{
			MrBayesPrint ("%s   WARNING: Report format for ratemult (parameter %d) varies across partitions.\n", spacer);
			MrBayesPrint ("%s      MrBayes will use the format for the first partition, which is %s.\n", spacer, modelParams[i].ratemultFormat);
			}
		}
	   
	/* probably a good idea to clean up link table here */
	paramCount = 0;
	for (j=0; j<NUM_LINKED; j++)
		{
		for (i=0; i<MAX_NUM_DIVS; i++)
			check[i] = NO;
		for (i=0; i<numCurrentDivisions; i++)
			{
			if (check[i] == NO && activeParams[j][i] > 0)
				{
				m = activeParams[j][i];
				paramCount++;
				for (k=i; k<numCurrentDivisions; k++)
					{
					if (check[k] == NO && activeParams[j][k] == m)
						{
						activeParams[j][k] = paramCount;
						check[k] = YES;
						}
					}
				}
			}
		}
				
	return (NO_ERROR);

}





int Unlink (void)

{
	
	return (NO_ERROR);

}
