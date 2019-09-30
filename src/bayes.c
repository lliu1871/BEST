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
#include <stdarg.h>
#include "mb.h"
#include "globals.h"
#include "bayes.h"
#include "command.h"
#include "mcmc.h"

#ifdef USE_READLINE
#include <readline/readline.h>
#include <readline/history.h>
static char **readline_completion(const char *, int, int);
#endif

#if defined (__MWERKS__)
#	if defined (MAC_VERSION)
#		include "SIOUX.h"
/* My codewarior version doesn't has this file, but compiles fine without Paul #		include "fonts.h" */
#	elif defined (WIN_VERSION)
#		include <windows.h>
#		include <wincon.h>
#		include <crtdbg.h>
#	endif
#endif
#if defined (WIN_VERSION)
#	undef NO_ERROR
#	undef ERROR
#	include <windows.h>
#	include <winbase.h>
#	undef NO_ERROR
#	undef ERROR
#	define NO_ERROR					0
#	define ERROR					1
#endif

/* local prototypes */
int  CommandLine (int argc, char **argv);
void PrintHeader (void);

/* global variables, declared in this file */
int			defTaxa;                     /* flag for whether number of taxa is known      */
int			defChars;                    /* flag for whether number of characters is known*/
int			defMatrix;                   /* flag for whether matrix is successfull read   */
int			defPartition;                /* flag for whether character partition is read  */
int			defConstraints;              /* flag for whether constraints on tree are read */
int			defPairs;                    /* flag for whether pairs are read               */
Doublet		doublet[16];                 /* holds information on states for doublets      */
int			fileNameChanged;			 /* has file name been changed ?                  */
long int	globalSeed;                  /* seed that is initialized at start up          */
int			nBitsInALong;                /* number of bits in a long                      */
int			readWord;					 /* should we read word next ?                    */
long int	runIDSeed;                   /* seed used only for determining run ID [stamp] */
long int	swapSeed;                    /* seed used only for determining which to swap  */
int         userLevel;                   /* user level                                    */

#			if defined (MPI_ENABLED)
int 		proc_id;                     /* process ID (0, 1, ..., num_procs-1)           */
int 		num_procs;                   /* number of active processors                   */
MrBFlt		myStateInfo[4];              /* likelihood/prior/heat vals of me              */
MrBFlt		partnerStateInfo[4];		 /* likelihood/prior/heat vals of partner         */
#			endif





int main (int argc, char *argv[])

{
	int i;

#	if defined (MPI_ENABLED)
	int		ierror;
#	endif

#	if defined (WIN_VERSION)
	HANDLE scbh;
	BOOL ok;
	DWORD lastError;
	COORD largestWindow;
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	int currBottom;
	char poltmp[256];

	scbh = GetStdHandle(STD_OUTPUT_HANDLE);
	GetConsoleScreenBufferInfo(scbh, &csbi);
	currBottom		= csbi.srWindow.Bottom;
	largestWindow   = GetLargestConsoleWindowSize(scbh);

	/* Allow for screen buffer 3000 lines long and 140 characters wide */
	csbi.dwSize.Y = 3000;
	csbi.dwSize.X = 140;

	SetConsoleScreenBufferSize(scbh, csbi.dwSize);
	/* Allow for maximum possible screen height */
	csbi.srWindow.Left		= 0; /* no change relative to current value */
	csbi.srWindow.Top		= 0; /* no change relative to current value */
	csbi.srWindow.Right		= 0; /* no change relative to current value */
	csbi.srWindow.Bottom	= largestWindow.Y - currBottom -10; /**/
	ok = SetConsoleWindowInfo(scbh, FALSE, &csbi.srWindow);
	if (ok == FALSE)
		{
		lastError = GetLastError();
		GetConsoleScreenBufferInfo(scbh, &csbi);
		sprintf(poltmp, "\nlastError = %d", lastError);
		printf(poltmp);
		}
#	endif
	/*mtrace();*/
	/* calculate the size of a long - used by bit manipulation functions */
	nBitsInALong = sizeof(long) * 8;
	if (nBitsInALong > 32) /* Do not use more than 32 bits until we    */
		nBitsInALong = 32; /* understand how 64-bit longs are handled. */

#	if defined (__MWERKS__) & defined (MAC_VERSION)
	/* Set up interface when using the Metrowerks compiler. This
	   should work for either Macintosh or Windows. */
	SIOUXSetTitle("\pMrBayes v3.1.2");
	SIOUXSettings.fontface         = 0;  /* plain=0; bold=1 */
	SIOUXSettings.setupmenus       = 0;
	SIOUXSettings.autocloseonquit  = 1;
	SIOUXSettings.asktosaveonclose = 0;
	SIOUXSettings.rows             = 60;
	SIOUXSettings.columns          = 90;
#	endif

#	if defined (MPI_ENABLED)
#		if defined (MAC_VERSION)
		printf ("                             Parallel version of\n\n");
#		else
		MrBayesPrint ("                             Parallel version of\n\n");
#		endif
	ierror = MPI_Init(&argc, &argv);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem initializing MPI\n", spacer);
		exit (1);
		}
	ierror = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem getting the number of processors\n", spacer);
		exit (1);
		}
	ierror = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem getting processors rank\n", spacer);
		exit (1);
		}
#	endif

#ifdef USE_READLINE
	rl_attempted_completion_function = readline_completion;
#endif
	/* Set up parameter table. */
	SetUpParms ();

	/* initialize seed using current time */
	GetTimeSeed ();

	/* Initialize the variables of the program. */
	InitializeMrBayes ();
	
	/* Print the nifty header. */
	PrintHeader ();
	
	/* Go to the command line, process any arguments passed to the program
	   and then wait for input. */
	i = CommandLine (argc, argv);

#	if defined (MPI_ENABLED)
	MPI_Finalize();
#	endif

	if (i == ERROR)
		return (1);
	else
		return (0);
	
}





int CommandLine (int argc, char **argv)

{

	int		i, message, nProcessedArgs;
	char	cmdStr[CMD_STRING_LENGTH];
#ifdef USE_READLINE
        char    *cmdStrP;
#endif
#			if defined (MPI_ENABLED)
	int		ierror;
#			endif

	for(i=0;i<CMD_STRING_LENGTH;i++) cmdStr[i]='\0';
	
	/* wait for user-input commands */
	nProcessedArgs = 1;	/* first argument is program name and needs not be processed */
	if (nProcessedArgs < argc)
		{
		mode = NONINTERACTIVE;	/* when a command is passed into the program, the default is to exit without listening to stdin */
		autoClose = YES;
		autoOverwrite = YES;
		noWarn = YES;
		quitOnError = YES;
		}
	for (;;)
		{
		if (nProcessedArgs < argc) 
			{
			/* we are here only if a command that has been passed
			   into the program remains to be processed */
			if (nProcessedArgs == 1 && (strcmp(argv[1],"-i") == 0 || strcmp(argv[1],"-I") == 0))
				{
				mode = INTERACTIVE;
				autoClose = NO;
				autoOverwrite = YES;
				noWarn = NO;
				quitOnError = NO;
				}
			else
				sprintf (cmdStr, "Execute %s", argv[nProcessedArgs]);
			nProcessedArgs++;
			}
		else
			{
			/* first check if we are in noninteractive mode and quit if so */
			if (mode == NONINTERACTIVE)
				{
				MrBayesPrint ("%s   Tasks completed, exiting program because mode is noninteractive\n", spacer);
				MrBayesPrint ("%s   To return control to the command line after completion of file processing, \n", spacer);
				MrBayesPrint ("%s   set mode to interactive with 'mb -i <filename>' (i is for interactive)\n", spacer);
				MrBayesPrint ("%s   or use 'set mode=interactive'\n\n", spacer);
				return (NO_ERROR);
				}
			/* normally, we simply wait at the prompt for a
			   user action */
#			if defined (MPI_ENABLED)
			if (proc_id == 0)
				{
	#ifdef USE_READLINE
		    cmdStrP = readline("Best > ");
			if(cmdStrP!=NULL) 
			        {
					strncpy(cmdStr,cmdStrP,CMD_STRING_LENGTH - 2);
					if (*cmdStrP) 
						add_history(cmdStrP);
					free (cmdStrP);
			        }
			else /* fall through to if (feof(stdin))..*/
     #else
				MrBayesPrint ("Best > ");
				fflush (stdin);
				if (fgets (cmdStr, CMD_STRING_LENGTH - 2, stdin) == NULL)
     #endif
					{
					if (feof(stdin))
						MrBayesPrint ("%s   End of File encountered on stdin; quitting\n", spacer);
					else
						MrBayesPrint ("%s   Could not read command from stdin; quitting\n", spacer);
					strcpy (cmdStr,"quit;\n");
					}
				}
			ierror = MPI_Bcast (&cmdStr, CMD_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
			if (ierror != MPI_SUCCESS)
				{
				MrBayesPrint ("%s   Problem broadcasting command string\n", spacer);
				}
#			else
	#ifdef USE_READLINE
		    cmdStrP = readline("Best > ");
			if(cmdStrP!=NULL) 
			        {
					strncpy(cmdStr,cmdStrP,CMD_STRING_LENGTH - 2);
					if (*cmdStrP) 
						add_history(cmdStrP);
					free (cmdStrP);
			        }
			else /* fall through to if (feof(stdin))..*/
	#else
			MrBayesPrint ("Best > ");
			fflush (stdin);
			if (fgets (cmdStr, CMD_STRING_LENGTH - 2, stdin) == NULL)
	#endif
				{
				if (feof(stdin))
					MrBayesPrint ("%s   End of File encountered on stdin; quitting\n", spacer);
				else
					MrBayesPrint ("%s   Could not read command from stdin; quitting\n", spacer);
				strcpy (cmdStr,"quit;\n");
				}
#			endif
			}
		i = 0;
		while (cmdStr[i] != '\0' && cmdStr[i] != '\n')
			i++;
		cmdStr[i++] = ';';
		cmdStr[i] = '\0';
		MrBayesPrint ("\n");
		if (cmdStr[0] != ';')
			{
			/* check that all characters in the string are valid */
			if (CheckStringValidity (cmdStr) == ERROR)
				{
				MrBayesPrint ("   Unknown character in command string\n\n");
				}
			else
				{
				expecting = Expecting(COMMAND);
				message = ParseCommand (cmdStr);

				if (message == NO_ERROR_QUIT)
					return (NO_ERROR);

				if (message == ERROR && quitOnError == YES)
					{
					MrBayesPrint ("%s   Will exit with signal 1 (error) because quitonerror is set to yes\n", spacer);
					MrBayesPrint ("%s   If you want control to be returned to the command line on error,\n", spacer);
					MrBayesPrint ("%s   use 'mb -i <filename>' (i is for interactive) or use 'set quitonerror=no'\n\n", spacer);
					return (ERROR);
					}			

#				if defined (MPI_ENABLED)
				ierror = MPI_Barrier (MPI_COMM_WORLD);
				if (ierror != MPI_SUCCESS)
					{
					MrBayesPrint ("%s   Problem at command barrier\n", spacer);
					}
#				endif

				MrBayesPrint ("\n");
				}
			}
		}

}

#ifdef USE_READLINE

extern char *command_generator(const char *text, int state);

char **readline_completion(const char *text, int start, int stop) {
	char **matches = (char **) NULL;
	
	if(start == 0)
    	matches = rl_completion_matches (text, command_generator);

  return (matches);	
}
#endif

void GetTimeSeed (void)

{

	time_t		curTime;

#	if defined (MPI_ENABLED)
	int			ierror;
	
	if (proc_id == 0)
		{
		curTime = time(NULL);
		globalSeed  = (long int)curTime;
		if (globalSeed < 0)
			globalSeed = -globalSeed;
		}
	ierror = MPI_Bcast(&globalSeed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem broadcasting seed\n", spacer);
		}
		
	if (proc_id == 0)
		{
		curTime = time(NULL);
		swapSeed  = (long int)curTime;
		if (swapSeed < 0)
			swapSeed = -swapSeed;
		}
	MPI_Bcast(&swapSeed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem broadcasting swap seed\n", spacer);
		}

	if (proc_id == 0)
		{
		curTime = time(NULL);
		runIDSeed  = (long int)curTime;
		if (runIDSeed < 0)
			runIDSeed = -runIDSeed;
		}
	ierror = MPI_Bcast(&runIDSeed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem broadcasting run ID seed\n", spacer);
		}		

#	else
	curTime = time(NULL);
	globalSeed  = (long int)curTime;
	if (globalSeed < 0)
		globalSeed = -globalSeed;
		
	curTime = time(NULL);
	swapSeed  = (long int)curTime;
	if (swapSeed < 0)
		swapSeed = -swapSeed;

	curTime = time(NULL);
	runIDSeed  = (long int)curTime;
	if (runIDSeed < 0)
		runIDSeed = -globalSeed;
		
#	endif
		
}





void InitializeMrBayes (void)

{
	/* this function initializes the program; only call it at the start of execution */
	
	int		i;

	userLevel            = STANDARD_USER;            /* default user level                            */

	readWord			 = NO;                       /* should we read a word next ?                  */
	fileNameChanged		 = NO;                       /* has the file name been changed ?              */
	echoMB                 = YES;                      /* flag used by Manual to control printing       */
	longIntegerSize      = sizeof(long int);         /* size of an long integer                       */

#	if defined (MPI_ENABLED)
	sprintf(manFileName, "commref_mb%sp.txt", VERSION_NUMBER);  /* name of command reference file      */
#	else
	sprintf(manFileName, "commref_mb%s.txt", VERSION_NUMBER); /* name of command reference file        */
#	endif

	for (i=0; i<NUM_ALLOCS; i++)                     /* set allocated memory to NO                    */
		memAllocs[i] = NO;				
	logToFile = NO;                                  /* should screen output be logged to a file      */
	strcpy(logFileName, "log.out");                  /* name of the log file                          */
	logFileFp = NULL;                                /* file pointer to log file                      */
	replaceLogFile = YES;                            /* should logfile be replace/appended to         */
	autoClose = NO;                                  /* set default autoclose                         */
	autoOverwrite = YES;                             /* set default autoOverwrite                     */
	noWarn = NO;                                     /* set default                             */
	quitOnError = NO;                                /* set default quitOnError                       */
	inferAncStates = NO;                             /* set default inferAncStates                    */
	inferSiteRates = NO;                             /* set default inferSiteRates                    */
	inferPosSel = NO;                                /* set default inferPosSel                       */
	inComment = NO;									 /* not in comment								  */
	numComments = 0;								 /* no comments encountered yet                   */
	mode = INTERACTIVE;								 /* set default mode							  */
	numOpenExeFiles = 0;							 /* no execute files open yet   				  */

	/* finally reset everything dependent on a matrix being defined */
	ReinitializeMrBayes ();
	
}





void PrintHeader (void)

{

#	if defined (MAC_VERSION)

#		if !defined (MPI_ENABLED)
		printf ("\n\n");
#		endif
		printf ("		     BEST V%S(Bayesian Estimation of Species Trees)\n\n", BESTVERSION_NUMBER);
		printf ("                                 Liang Liu\n\n" ); 
		printf ("                   Department of Organismic and Evolutionary Biology\n" );
		printf ("                                 Harvard University\n" );
		printf ("                                lliu@oeb.harvard.edu\n\n\n\n" );   
		printf ("                            based on MrBayes v%s\n\n", VERSION_NUMBER);
		printf ("                      (Bayesian Analysis of Phylogeny)\n\n");
#		if defined (MPI_ENABLED)
		printf ("                             (Parallel version)\n");
		printf ("                         (%d processors available)\n\n", num_procs);
#		endif
		srand ((unsigned int)time(NULL));
		if (rand() % 2)
			{
			printf ("                                     by\n\n");
			printf ("                  John P. Huelsenbeck and Fredrik Ronquist\n\n");
			printf ("                 Section of Ecology, Behavior and Evolution\n");
			printf ("                       Division of Biological Sciences\n");
			printf ("                     University of California, San Diego\n");
			printf ("                           johnh@biomail.ucsd.edu\n\n");
			printf ("                       School of Computational Science\n");
			printf ("                           Florida State University\n");
			printf ("                            ronquist@csit.fsu.edu \n\n");
			}
		else
			{
			printf ("                                     by\n\n");
			printf ("                  Fredrik Ronquist and John P. Huelsenbeck\n\n");
			printf ("                       School of Computational Science\n");
			printf ("                           Florida State University\n");
			printf ("                            ronquist@csit.fsu.edu \n\n");
			printf ("                 Section of Ecology, Behavior and Evolution\n");
			printf ("                       Division of Biological Sciences\n");
			printf ("                     University of California, San Diego\n");
			printf ("                           johnh@biomail.ucsd.edu\n\n");
			}
		printf ("              Distributed under the GNU General Public License\n\n");

#		if defined (MPI_ENABLED)
		if (proc_id == 0)
			{
			printf ("                               Node number %d\n\n", proc_id + 1);
			printf ("               Type \"help\" or \"help <command>\" for information\n");
			printf ("                     on the commands that are available.\n\n\n");
			}
		else
			{
			printf ("                           Remote node number %d\n\n", proc_id + 1);
			printf ("                    All information is printed to node 1\n\n\n");
			}
#		else
		printf ("               Type \"help\" or \"help <command>\" for information\n");
		printf ("                     on the commands that are available.\n\n\n");
#		endif

#	else       
               
#		if !defined (MPI_ENABLED)
		MrBayesPrint ("\n\n");
#		endif
		MrBayesPrint ("                   BEST v%s(Bayesian Estimation of Species Trees)\n\n",BESTVERSION_NUMBER);
                MrBayesPrint ("                                 Liang Liu\n\n");
		MrBayesPrint ("                   Department of Organismic and Evolutionary Biology\n" );              
                MrBayesPrint ("                                 Harvard University\n" );  
                MrBayesPrint ("                                lliu@oeb.harvard.edu\n\n\n\n" ); 
		MrBayesPrint ("                            based on MrBayes v%s\n\n", VERSION_NUMBER);
		MrBayesPrint ("                      (Bayesian Analysis of Phylogeny)\n\n");
#		if defined (MPI_ENABLED)
		MrBayesPrint ("                             (Parallel version)\n");
		MrBayesPrint ("                         (%d processors available)\n\n", num_procs);
#		endif
		srand((unsigned int) time (NULL));
		if (rand() % 2)
			{
			MrBayesPrint ("                                     by\n\n");
			MrBayesPrint ("                  John P. Huelsenbeck and Fredrik Ronquist\n\n");
			MrBayesPrint ("                 Section of Ecology, Behavior and Evolution\n");
			MrBayesPrint ("                       Division of Biological Sciences\n");
			MrBayesPrint ("                     University of California, San Diego\n");
			MrBayesPrint ("                           johnh@biomail.ucsd.edu\n\n");
			MrBayesPrint ("                       School of Computational Science\n");
			MrBayesPrint ("                           Florida State University\n");
			MrBayesPrint ("                            ronquist@csit.fsu.edu \n\n");
			}
		else
			{
			MrBayesPrint ("                                     by\n\n");
			MrBayesPrint ("                  Fredrik Ronquist and John P. Huelsenbeck\n\n");
			MrBayesPrint ("                       School of Computational Science\n");
			MrBayesPrint ("                           Florida State University\n");
			MrBayesPrint ("                            ronquist@csit.fsu.edu \n\n");
			MrBayesPrint ("                 Section of Ecology, Behavior and Evolution\n");
			MrBayesPrint ("                       Division of Biological Sciences\n");
			MrBayesPrint ("                     University of California, San Diego\n");
			MrBayesPrint ("                           johnh@biomail.ucsd.edu\n\n");
			}
		MrBayesPrint ("              Distributed under the GNU General Public License\n\n");
		MrBayesPrint ("               Type \"help\" or \"help <command>\" for information\n");
		MrBayesPrint ("                     on the commands that are available.\n\n\n");	
#	endif
	
}





int ReinitializeMrBayes (void)

{

	/* this function resets everything dependent on a matrix */
	
	int				i;
	
	/* reset all matrix flags */
	defTaxa              = NO;                       /* flag for whether number of taxa is known      */
	defChars             = NO;                       /* flag for whether number of characters is known*/
	defMatrix            = NO;                       /* flag for whether matrix is successfull read   */
	matrixHasPoly		 = NO;                       /* flag for whether matrix has polymorphisms     */
    isInAmbig			 = NO;						 /* flag for whether the parser is within ()      */
	isInPoly			 = NO;						 /* flag for whether the parser is within {}      */
	defPartition         = NO;                       /* flag for whether character partition is read  */
	defConstraints       = NO;                       /* flag for whether constraints on tree are read */
	defPairs             = NO;                       /* flag indicating whether pairs have been defnd */
	numDefinedPartitions = 0;                        /* number of defined partitions                  */
	partitionNum         = 1;                        /* partition number currently enforced           */
	numCurrentDivisions  = 0;                        /* number of partitions of data                  */
	numCharSets          = 0;          			     /* holds number of character sets                */
	numPartitions		 = 1;          			     /* holds number of partitions                    */
	numDefinedConstraints= 0;          			     /* holds number of defined constraints           */
	isTranslateDef       = NO;          			 /* is translate block defined?                   */
	outGroupNum			 = 0;            			 /* default outgroup                              */
	isMixed			     = NO;          			 /* are data mixed ?                              */
	isUserTreeDefined    = NO;          			 /* is user tree defined?                         */
	dataType             = NONE;          			 /* holds datatype                                */
	matchId = gapId = missingId = 'z';

	if (InitializeLinks () == ERROR)
		{
		MrBayesPrint ("%s   Problem initializing link table\n", spacer);
		return (ERROR);
		}

	SetModelDefaults ();

	strcpy (spacer, "");                             /* holds blanks for indentation                  */

	/* chain parameters */
	chainParams.numGen = 1000000;                    /* number of MCMC cycles                         */
	chainParams.sampleFreq = 100;                    /* frequency to sample chain                     */
	chainParams.printFreq = 100;                     /* frequency to print chain                      */
	chainParams.swapFreq = 1;                        /* frequency of attempting swap of states        */
	chainParams.numSwaps = 1;                        /* number of swaps to try each time              */
	chainParams.mcmcDiagn = YES;                     /* write MCMC diagnostics to file ?              */
	chainParams.diagnFreq = 1000;                    /* diagnostics frequency                         */
	chainParams.minPartFreq = 0.10;                  /* min partition frequency for diagnostics       */
	chainParams.allChains = NO;                      /* calculate diagnostics for all chains ?        */
	chainParams.allComps = NO;                       /* do not calc diagn for all run comparisons     */
	chainParams.relativeBurnin = YES;                /* use relative burnin?                          */
	chainParams.burninFraction = 0.25;				 /* default burnin fraction                       */
	chainParams.stopRule = NO;						 /* should stopping rule be used?                 */
	chainParams.stopVal = 0.01;						 /* convergence diagnostic value to reach         */
	chainParams.numRuns = 2;                         /* number of runs                                */
	chainParams.numChains = 4;                       /* number of chains                              */
	chainParams.chainTemp = 0.2;                     /* chain temperature                             */
	chainParams.redirect = NO;                       /* should printf be to stdout                    */
	strcpy(chainParams.chainFileName, "temp.out");   /* chain file name for output                    */
	chainParams.chainBurnIn = 0;                     /* chain burn in length                          */
	chainParams.numStartPerts = 0;                   /* number of perturbations to starting tree      */
	strcpy(chainParams.chainStartTree, "Random");    /* starting tree for chain (random/user)         */
	chainParams.saveBrlens = YES;                    /* should branch lengths be saved                */
	chainParams.chainSeed = globalSeed;              /* random seed for chain                         */
	chainParams.weightScheme[0] = 0.0;               /* percent chars to decrease in weight           */
	chainParams.weightScheme[1] = 0.0;               /* percent chars to increase in weight           */
	chainParams.weightScheme[2] = 1.0;               /* weight increment                              */
	chainParams.calcPbf = NO;                        /* should we calculate the pseudo-BF?            */
	chainParams.pbfInitBurnin = 100000;              /* initial burnin for pseudo BF                  */
	chainParams.pbfSampleFreq = 10;                  /* sample frequency for pseudo BF                */
	chainParams.pbfSampleTime = 2000;                /* how many cycles to calcualate site prob.      */
	chainParams.pbfSampleBurnin = 2000;              /* burnin period for each site for pseudo BF     */
	chainParams.userDefinedTemps = NO;               /* should we use the users temperatures?         */
	for (i=0; i<MAX_CHAINS; i++)
		chainParams.userTemps[i] = 1.0;              /* user-defined chain temperatures               */
	chainParams.swapAdjacentOnly = NO;               /* swap only adjacent temperatures               */
	chainParams.printMax = 8;						 /* maximum number of chains to print to screen   */
	chainParams.printAll = YES;						 /* whether to print heated chains                */
	chainParams.runWithData = YES;				     /* whether to run with data                      */
	chainParams.orderTaxa = NO;                      /* should taxa be ordered in output trees?       */
	strcpy (chainParams.mpiPath, "");      		 /* path to temporary files for mpi version       */

	/* sumt parameters */
	strcpy(sumtParams.sumtFileName, "temp.t");       /* input name for sumt command                   */
	sumtParams.sumtBurnIn = 0;                       /* burnin for sumt command                       */
	strcpy(sumtParams.sumtConType, "Halfcompat");    /* type of consensus tree output                 */
	sumtParams.calcTrprobs = YES;                    /* should individual tree probs be calculated    */
	sumtParams.showSumtTrees = NO;                   /* should the individual tree probs be shown     */
	sumtParams.freqDisplay = 0.05;                   /* threshold for printing partitions to screen   */
	sumtParams.printBrlensToFile = NO;               /* should brlens be printed to file              */
	sumtParams.brlensFreqDisplay = 0.20;             /* threshold for printing brlens to file         */
	sumtParams.numTrees = 1;                         /* number of trees to summarize                  */
	sumtParams.numRuns = 2;                          /* number of analyses to summarize               */
	sumtParams.orderTaxa = YES;						 /* order taxa in trees ?                         */

	/* sump parameters */
	strcpy(sumpParams.sumpFileName, "temp.p");       /* input name for sump command                   */
	sumpParams.sumpBurnIn = 0;                       /* burnin for sump command                       */
	sumpParams.margLike = YES;					     /* print marginal likelihood                     */
	sumpParams.plot = YES;                           /* print likelihood plot                         */
	sumpParams.table = YES;                          /* print parameter table                         */
	sumpParams.printToFile = NO;                     /* do not print output to file                   */
	strcpy (sumpParams.sumpOutfile, "temp.p.stat");  /* output name for sump command                  */
	sumpParams.numRuns = 2;                          /* number of analyses to summarize               */

	/* comparetree parameters */
	strcpy(comptreeParams.comptFileName1, "temp.t"); /* input name for comparetree command            */
	strcpy(comptreeParams.comptFileName2, "temp.t"); /* input name for comparetree command            */
	strcpy(comptreeParams.comptOutfile, "temp.comp");/* input name for comparetree command            */
	comptreeParams.comptBurnIn = 0;                  /* burnin for comparetree command                */

	/* plot parameters */
	strcpy(plotParams.plotFileName, "temp.p");       /* input name for plot command                   */
	strcpy(plotParams.parameter, "lnL");             /* plotted parameter plot command                */
	strcpy(plotParams.match, "Perfect");             /* matching for plot command                     */
	plotParams.plotBurnIn = 0;                       /* burnin for plot command                       */
	
	/* set the proposal information */
	SetUpMoveTypes ();

	/* set up doublet information */
	doublet[ 0].first  = 1;   doublet[ 0].second = 1;
	doublet[ 1].first  = 1;   doublet[ 1].second = 2;
	doublet[ 2].first  = 1;   doublet[ 2].second = 4;
	doublet[ 3].first  = 1;   doublet[ 3].second = 8;
	doublet[ 4].first  = 2;   doublet[ 4].second = 1;
	doublet[ 5].first  = 2;   doublet[ 5].second = 2;
	doublet[ 6].first  = 2;   doublet[ 6].second = 4;
	doublet[ 7].first  = 2;   doublet[ 7].second = 8;
	doublet[ 8].first  = 4;   doublet[ 8].second = 1;
	doublet[ 9].first  = 4;   doublet[ 9].second = 2;
	doublet[10].first  = 4;   doublet[10].second = 4;
	doublet[11].first  = 4;   doublet[11].second = 8;
	doublet[12].first  = 8;   doublet[12].second = 1;
	doublet[13].first  = 8;   doublet[13].second = 2;
	doublet[14].first  = 8;   doublet[14].second = 4;
	doublet[15].first  = 8;   doublet[15].second = 8;

	return (NO_ERROR);

}





int DoQuit (void)

{

	int			i;
	char		tempName[100];
		
	/* free information for matrix */
	FreeMatrix ();
	
	SafeFclose (&logFileFp);

	/* check to see if any memory has not been freed */
	for (i=0; i<NUM_ALLOCS; i++)
		{
		if (memAllocs[i] == YES)
			{
			MrBayesPrint ("   WARNING: Memory (%d) has not been freed\n", i);
			if (mode == INTERACTIVE && quitOnError == NO)
				{
				MrBayesPrint ("%s   Hit return key to continue  ", spacer);
				fflush (stdin);
				fgets (tempName, 100, stdin);
				}
			}
		}

	MrBayesPrint ("   Quitting program\n\n");

	/* If we quit while reading a mrbayes block, then we need to make certain
	   that we return a NO_ERROR_QUIT so we can break out of DoExecute cleanly,
	   and dealloc "s" there. */
	if (inMrbayesBlock == YES)
		{
		inMrbayesBlock = NO;
		return (NO_ERROR_QUIT);
		}
		
	return (NO_ERROR);

}





void SetModelDefaults (void)

{

	int			i, j;

	/* model parameters */
	for (j=0; j<MAX_NUM_DIVS; j++)
		{
		modelParams[j].dataType = -1;                     /* data type for partition                      */
		strcpy(modelParams[j].nucModel, "4by4");          /* nucleotide model used                        */
		strcpy(modelParams[j].nst, "1");                  /* number of substitution types                 */
		strcpy(modelParams[j].aaModelPr, "Fixed");        /* amino acid model prior used                  */
		strcpy(modelParams[j].aaModel, "Poisson");        /* amino acid model used                        */
		for (i=0; i<10; i++)
			modelParams[j].aaModelPrProbs[i] = 0.0;
		strcpy(modelParams[j].parsModel, "No");           /* use the (so-called) parsimony model          */
		strcpy(modelParams[j].geneticCode, "Universal");  /* genetic code used                            */
		SetCode (j);
		strcpy(modelParams[j].ploidy, "Diploid");         /* ploidy level                                 */
		strcpy(modelParams[j].coding, "All");             /* type of patterns encoded                     */
		strcpy(modelParams[j].omegaVar, "Equal");         /* type of omega variation model                */
		strcpy(modelParams[j].ratesModel, "Equal");       /* rates across sites model                     */
		modelParams[j].numGammaCats = 4;                  /* number of categories for gamma approximation */
		modelParams[j].numBetaCats = 5;                   /* number of categories for beta approximation  */
		strcpy(modelParams[j].covarionModel, "No");       /* use covarion model? (yes/no)                 */
		strcpy(modelParams[j].augmentData, "No");         /* should data be augmented                     */
		strcpy(modelParams[j].tRatioPr, "Beta");          /* prior for ti/tv rate ratio                   */
		modelParams[j].tRatioFix = 1.0;
		modelParams[j].tRatioDir[0] = 1.0;
		modelParams[j].tRatioDir[1] = 1.0;
		strcpy(modelParams[j].revMatPr, "Dirichlet");     /* prior for GTR model (nucleotides)            */
		for (i=0; i<6; i++)
			{
			modelParams[j].revMatFix[i] = 1.0;
			modelParams[j].revMatDir[i] = 1.0;
			}
		strcpy(modelParams[j].aaRevMatPr, "Dirichlet");     /* prior for GTR model (proteins)             */
		for (i=0; i<190; i++)
			{
			modelParams[j].aaRevMatFix[i] = 1.0;
			modelParams[j].aaRevMatDir[i] = 1.0;
			}
		strcpy(modelParams[j].omegaPr, "Dirichlet");      /* prior for omega                              */
		modelParams[j].omegaFix = 1.0;
		modelParams[j].omegaDir[0] = 1.0;
		modelParams[j].omegaDir[1] = 1.0;
		strcpy(modelParams[j].ny98omega1pr, "Beta");      /* prior for class 1 omega (Ny98 model)         */
		modelParams[j].ny98omega1Fixed = 0.1;
		modelParams[j].ny98omega1Beta[0] = 1.0;
		modelParams[j].ny98omega1Beta[1] = 1.0;
		strcpy(modelParams[j].ny98omega3pr, "Exponential");/* prior for class 3 omega (Ny98 model)        */
		modelParams[j].ny98omega3Fixed = 2.0;
		modelParams[j].ny98omega3Uni[0] = 1.0;
		modelParams[j].ny98omega3Uni[1] = 50.0;
		modelParams[j].ny98omega3Exp = 1.0;
		strcpy(modelParams[j].m3omegapr, "Exponential");  /* prior for all three omegas (M3 model)        */
		modelParams[j].m3omegaFixed[0] = 0.1;
		modelParams[j].m3omegaFixed[1] = 1.0;
		modelParams[j].m3omegaFixed[2] = 2.0;
		strcpy(modelParams[j].m10betapr, "Uniform");      /* prior for omega variation (M10 model)        */
		strcpy(modelParams[j].m10gammapr, "Uniform");
		modelParams[j].m10betaUni[0] = 0.0;
		modelParams[j].m10betaUni[1] = 20.0;
		modelParams[j].m10betaExp = 1.0;
		modelParams[j].m10betaFix[0] = 1.0;
		modelParams[j].m10betaFix[1] = 1.0;
		modelParams[j].m10gammaUni[0] = 0.0;
		modelParams[j].m10gammaUni[1] = 20.0;
		modelParams[j].m10gammaExp = 1.0;
		modelParams[j].m10gammaFix[0] = 1.0;
		modelParams[j].m10gammaFix[1] = 1.0;
		modelParams[j].numM10GammaCats = 4;
		modelParams[j].numM10BetaCats = 4;
		strcpy(modelParams[j].codonCatFreqPr, "Dirichlet");/* prior for selection cat frequencies         */
		modelParams[j].codonCatFreqFix[0] = 1.0/3.0;
		modelParams[j].codonCatFreqFix[1] = 1.0/3.0;
		modelParams[j].codonCatFreqFix[2] = 1.0/3.0;
		modelParams[j].codonCatDir[0] = 1.0;
		modelParams[j].codonCatDir[1] = 1.0;
		modelParams[j].codonCatDir[2] = 1.0;
		strcpy(modelParams[j].stateFreqPr, "Dirichlet");  /* prior for character state frequencies        */
		strcpy(modelParams[j].stateFreqsFixType, "Equal");
		for (i=0; i<200; i++)
			{
			modelParams[j].stateFreqsFix[i] = 0.0;   
			modelParams[j].stateFreqsDir[i] = 1.0;
			}    
		modelParams[j].numDirParams = 0;
		strcpy(modelParams[j].shapePr, "Uniform");        /* prior for gamma shape parameter              */
		modelParams[j].shapeFix = 0.5;
		modelParams[j].shapeUni[0] = MIN_SHAPE_PARAM;
		modelParams[j].shapeUni[1] = MAX_SHAPE_PARAM;
		modelParams[j].shapeExp = 2.0;
		strcpy(modelParams[j].pInvarPr, "Uniform");       /* prior for proportion of invariable sites     */
		modelParams[j].pInvarFix = 0.1;
		modelParams[j].pInvarUni[0] = 0.0;
		modelParams[j].pInvarUni[1] = 1.0;
		strcpy(modelParams[j].adGammaCorPr, "Uniform");   /* prior for correlation param of adGamma model */
		modelParams[j].corrFix = 0.0;
		modelParams[j].corrUni[0] = -1.0;
		modelParams[j].corrUni[1] = 1.0;
		strcpy(modelParams[j].covSwitchPr, "Uniform");    /* prior for switching rates of covarion model  */
		modelParams[j].covswitchFix[0] = 1.0;
		modelParams[j].covswitchFix[1] = 1.0;
		modelParams[j].covswitchUni[0] = 0.0;
		modelParams[j].covswitchUni[1] = 100.0;
		modelParams[j].covswitchExp = 1.0;
		strcpy(modelParams[j].symPiPr, "Fixed");           /* prior for pi when unidentifiable states used */
		modelParams[j].symBetaFix = -1.0;
		modelParams[j].symBetaUni[0] = 0.0;
		modelParams[j].symBetaUni[1] = 20.0;
		modelParams[j].symBetaExp = 2;
		strcpy(modelParams[j].ratePr, "Fixed");           /* prior on rate for a partition                */
		modelParams[j].ratePrDir = 1.0;
		strcpy(modelParams[j].topologyPr, "Uniform");     /* prior for tree topology                      */
		for (i=0; i<30; i++)
			modelParams[j].activeConstraints[i] = NO;     /* which constraints are active                 */
		strcpy(modelParams[j].brlensPr, "Unconstrained"); /* prior on branch lengths                      */
		modelParams[j].brlensUni[0] = BRLENS_MIN;
		modelParams[j].brlensUni[1] = 10.0;
		modelParams[j].brlensExp = 10.0;
		strcpy(modelParams[j].unconstrainedPr, "Exponential");/* prior on branches if unconstrained           */
		strcpy(modelParams[j].clockPr, "Uniform");        /* prior on branch lengths if clock enforced    */
		strcpy(modelParams[j].speciationPr, "Uniform");   /* prior on speciation rate                     */
		modelParams[j].speciationFix = 1.0;
		modelParams[j].speciationUni[0] = 0.0;
		modelParams[j].speciationUni[1] = 10.0;
		modelParams[j].speciationExp = 1.0;
		strcpy(modelParams[j].extinctionPr, "Uniform");   /* prior on extinction rate                     */
		modelParams[j].extinctionFix = 1.0;
		modelParams[j].extinctionUni[0] = 0.0;
		modelParams[j].extinctionUni[1] = 10.0;
		modelParams[j].extinctionExp = 1.0;
		strcpy(modelParams[j].treeHeightPr, "Exponential");/* prior on tree height                        */
		modelParams[j].treeHeightGamma[0] = 1.0;
		modelParams[j].treeHeightGamma[1] = 1.0;
		modelParams[j].treeHeightExp = 1.0;
		modelParams[j].sampleProb = 1.0;                  /* taxon sampling fraction                      */
		strcpy(modelParams[j].thetaPr, "Uniform");        /* prior on coalescence prior                   */
		modelParams[j].thetaFix = 1.0;
		modelParams[j].thetaUni[0] = 0.0;
		modelParams[j].thetaUni[1] = 10.0;
		modelParams[j].thetaExp = 1.0;
		modelParams[j].thetaGamma[0] = 3.0;   /*prior for joint gene tree*/
		modelParams[j].thetaGamma[1] = 0.003;
	
		strcpy(modelParams[j].GeneMuPr, "Fixed");        /* prior on coalescence prior                   */
		modelParams[j].GeneMuFix = 1.0;

		strcpy(modelParams[j].growthPr, "Fixed");        /* prior on coalescence growth rate prior      */
		modelParams[j].growthFix = 0.0;
		modelParams[j].growthUni[0] = 0.0;
		modelParams[j].growthUni[1] = 100.0;
		modelParams[j].growthExp = 1.0;
		modelParams[j].growthNorm[0] = 0.0;
		modelParams[j].growthNorm[1] = 1.0;
		strcpy(modelParams[j].brownCorPr, "Fixed");       /* prior on correlation of brownian model       */
		modelParams[j].brownCorrFix = 0.0;
		modelParams[j].brownCorrUni[0] = -1.0;
		modelParams[j].brownCorrUni[1] = 1.0;
		strcpy(modelParams[j].brownScalesPr, "Gammamean");/* prior on scales of brownian model            */
		modelParams[j].brownScalesFix = 10.0;
		modelParams[j].brownScalesUni[0] = 0.0;
		modelParams[j].brownScalesUni[1] = 100.0;
		modelParams[j].brownScalesGamma[0] = 1.0;
		modelParams[j].brownScalesGamma[1] = 10.0;
		modelParams[j].brownScalesGammaMean = 10.0;
		strcpy (modelParams[j].calWaitPr, "Uniform") ;	/* prior on clock rate						 */
		modelParams[j].calWaitExp = 1.0 / 10.0;
		modelParams[j].calWaitFix = 1.0;
		modelParams[j].calWaitUni[0] = 0.0;
		modelParams[j].calWaitUni[1] = 10000000.0;

		strcpy(modelParams[j].tratioFormat, "Ratio");	    /* Default format for tratio				  */
		strcpy(modelParams[j].revmatFormat, "Dirichlet");	/* Default format for revmat				  */
		strcpy(modelParams[j].ratemultFormat, "Scaled");	/* Default format for ratemult				  */
		strcpy(modelParams[j].inferAncStates, "No");		/* Do not infer ancestral states			  */
		strcpy(modelParams[j].inferPosSel, "No");			/* Do not infer positive selection			  */
		strcpy(modelParams[j].inferSiteRates, "No");		/* Do not infer site rates					  */

		}

}


void *SafeMalloc(size_t s) {
        void *ptr = malloc(s);
        if(ptr==NULL)
                return NULL;
        return memset(ptr,0,s);
}

int SafeFclose(FILE **fp) {
	int retval=-1;
	if( fp!=NULL && (*fp)!=NULL ) 
		retval=fclose(*fp);
	*fp = NULL;
	return retval;	
}

void SetCode (int part)

{

	int			i, s, s1, s2, s3, ns;

	modelParams[part].codon[ 0] = 12; /* AAA Lys */
	modelParams[part].codon[ 1] =  3; /* AAC Asn */
	modelParams[part].codon[ 2] = 12; /* AAG Lys */
	modelParams[part].codon[ 3] =  3; /* AAT Asn */
	modelParams[part].codon[ 4] = 17; /* ACA Thr */
	modelParams[part].codon[ 5] = 17; /* ACC Thr */
	modelParams[part].codon[ 6] = 17; /* ACG Thr */
	modelParams[part].codon[ 7] = 17; /* ACT Thr */
	modelParams[part].codon[ 8] =  2; /* AGA Arg */
	modelParams[part].codon[ 9] = 16; /* AGC Ser */
	modelParams[part].codon[10] =  2; /* AGG Arg */
	modelParams[part].codon[11] = 16; /* AGT Ser */
	modelParams[part].codon[12] = 10; /* ATA Ile */
	modelParams[part].codon[13] = 10; /* ATC Ile */
	modelParams[part].codon[14] = 13; /* ATG Met */
	modelParams[part].codon[15] = 10; /* ATT Ile */
	modelParams[part].codon[16] =  6; /* CAA Gln */
	modelParams[part].codon[17] =  9; /* CAC His */
	modelParams[part].codon[18] =  6; /* CAG Gln */
	modelParams[part].codon[19] =  9; /* CAT His */
	modelParams[part].codon[20] = 15; /* CCA Pro */
	modelParams[part].codon[21] = 15; /* CCC Pro */
	modelParams[part].codon[22] = 15; /* CCG Pro */
	modelParams[part].codon[23] = 15; /* CCT Pro */
	modelParams[part].codon[24] =  2; /* CGA Arg */
	modelParams[part].codon[25] =  2; /* CGC Arg */
	modelParams[part].codon[26] =  2; /* CGG Arg */
	modelParams[part].codon[27] =  2; /* CGT Arg */
	modelParams[part].codon[28] = 11; /* CTA Leu */
	modelParams[part].codon[29] = 11; /* CTC Leu */
	modelParams[part].codon[30] = 11; /* CTG Leu */
	modelParams[part].codon[31] = 11; /* CTT Leu */
	modelParams[part].codon[32] =  7; /* GAA Glu */
	modelParams[part].codon[33] =  4; /* GAC Asp */
	modelParams[part].codon[34] =  7; /* GAG Glu */
	modelParams[part].codon[35] =  4; /* GAT Asp */
	modelParams[part].codon[36] =  1; /* GCA Ala */
	modelParams[part].codon[37] =  1; /* GCC Ala */
	modelParams[part].codon[38] =  1; /* GCG Ala */
	modelParams[part].codon[39] =  1; /* GCT Ala */
	modelParams[part].codon[40] =  8; /* GGA Gly */
	modelParams[part].codon[41] =  8; /* GGC Gly */
	modelParams[part].codon[42] =  8; /* GGG Gly */
	modelParams[part].codon[43] =  8; /* GGT Gly */
	modelParams[part].codon[44] = 20; /* GTA Val */
	modelParams[part].codon[45] = 20; /* GTC Val */
	modelParams[part].codon[46] = 20; /* GTG Val */
	modelParams[part].codon[47] = 20; /* GTT Val */
	modelParams[part].codon[48] = 21; /* TAA Stop*/
	modelParams[part].codon[49] = 19; /* TAC Tyr */
	modelParams[part].codon[50] = 21; /* TAG Stop*/
	modelParams[part].codon[51] = 19; /* TAT Tyr */
	modelParams[part].codon[52] = 16; /* TCA Ser */
	modelParams[part].codon[53] = 16; /* TCC Ser */
	modelParams[part].codon[54] = 16; /* TCG Ser */
	modelParams[part].codon[55] = 16; /* TCT Ser */
	modelParams[part].codon[56] = 21; /* TGA Stop*/
	modelParams[part].codon[57] =  5; /* TGC Cys */
	modelParams[part].codon[58] = 18; /* TGG Trp */
	modelParams[part].codon[59] =  5; /* TGT Cys */
	modelParams[part].codon[60] = 11; /* TTA Leu */
	modelParams[part].codon[61] = 14; /* TTC Phe */
	modelParams[part].codon[62] = 11; /* TTG Leu */
	modelParams[part].codon[63] = 14; /* TTT Phe */
	
	if (!strcmp(modelParams[part].geneticCode, "Vertmt"))
		{
		/* UGA: Ter -> Trp
		   AUA: Ile -> Met
		   AGA: Arg -> Ter 
		   AGG: Arg -> Ter */
		modelParams[part].codon[ 8] = 21; /* AGA Stop */ 
		modelParams[part].codon[10] = 21; /* AGG Stop */
		modelParams[part].codon[12] = 13; /* ATA Met  */
		modelParams[part].codon[56] = 18; /* TGA Trp  */
		}
	else if (!strcmp(modelParams[part].geneticCode, "Mycoplasma"))
		{
		/* UGA: Ter -> Trp */
		modelParams[part].codon[56] = 18; /* TGA Trp */
		}
	else if (!strcmp(modelParams[part].geneticCode, "Yeast"))
		{
		/* UGA: Ter -> Trp
		   AUA: Ile -> Met
		   CUA: Leu -> Thr
		   CUC: Leu -> Thr
		   CUG: Leu -> Thr
		   CUU: Leu -> Thr */
		modelParams[part].codon[12] = 13; /* ATA Met */
		modelParams[part].codon[28] = 17; /* CTA Thr */
		modelParams[part].codon[29] = 17; /* CTC Thr */
		modelParams[part].codon[30] = 17; /* CTG Thr */
		modelParams[part].codon[31] = 17; /* CTT Thr */
		modelParams[part].codon[56] = 18; /* TGA Trp */
		}
	else if (!strcmp(modelParams[part].geneticCode, "Ciliates"))
		{
		/* UAA: Ter -> Gln
		   UAG: Ter -> Gln */
		modelParams[part].codon[48] =  6; /* TAA Gln */
		modelParams[part].codon[50] =  6; /* TAG Gln */
		}
	else if (!strcmp(modelParams[part].geneticCode, "Metmt"))
		{
		/* UGA: Ter -> Trp
		   AUA: Ile -> Met
		   AGA: Arg -> Ser
		   AGG: Arg -> Ser */
		modelParams[part].codon[ 8] = 16; /* AGA Ser */ 
		modelParams[part].codon[10] = 16; /* AGG Ser */
		modelParams[part].codon[12] = 13; /* ATA Met */
		modelParams[part].codon[56] = 18; /* TGA Trp */
		}
	else
		{

		}
	
	ns = 0;
	for (i=0; i<64; i++)
		{
		if (modelParams[part].codon[i] != 21)
			ns++;
		}
	/* printf ("ns = %d\n", ns); */
	
	s = 0;
	for (s1=0; s1<4; s1++)
		{
		for (s2=0; s2<4; s2++)
			{
			for (s3=0; s3<4; s3++)
				{
				if (modelParams[part].codon[s1*16 + s2*4 + s3] != 21)
					{
					modelParams[part].codonNucs[s][0] = s1;
					modelParams[part].codonNucs[s][1] = s2;
					modelParams[part].codonNucs[s][2] = s3;
					modelParams[part].codonAAs[s] = modelParams[part].codon[s1*16 + s2*4 + s3];
					s++;
					}
				}
			}
		}
		
}





void MrBayesPrint (char *format, ...)

{
	va_list ptr;

#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
		if (echoMB == YES)
			{
			va_start (ptr, format);
			vprintf (format, ptr);
			va_end(ptr);
			fflush (stdout);
			}
		if (logToFile == YES)
			{
			if (logFileFp == NULL)
				printf ("%s   Could not print log output to file\n", spacer);
			else
				{
				va_start (ptr, format);
				vfprintf (logFileFp, format, ptr);
				va_end(ptr);
				fflush (logFileFp);
				}
			}
		}
#	else
	if (chainParams.redirect == NO)
		{
		if (echoMB == YES)
			{
			va_start (ptr, format);
			vprintf (format, ptr);
			va_end(ptr);
			fflush (stdout);
			}
		if (logToFile == YES)
			{
			if (logFileFp == NULL)
				{
				printf ("%s   Could not print log output to file\n", spacer);
				logToFile = NO;
				}
			else
				{
				va_start (ptr, format);
				vfprintf (logFileFp, format, ptr);
				va_end(ptr);
				fflush (logFileFp);
				}
			}
		}
#	endif
}




void MrBayesPrintf (FILE *f, char *format, ...)

{
	va_list                 ptr;

#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
		va_start (ptr, format);
		vfprintf (f, format, ptr);
		va_end(ptr);
		fflush(f);
		}
#	else
	va_start (ptr, format);
	vfprintf (f, format, ptr);
	va_end(ptr);
	fflush(f);
#	endif
}


