/* global variables */
extern int 			NumRuns;
extern SPTree			*gtree;
extern SPTree			sptree;
extern char		*printString;           /* string for printing to a file                */
extern size_t	printStringSize;       	/* length of printString                        */
extern double			poissonMean;
extern int              	jointGenePr;
extern int			spMupr,spMupralpha;
extern double			speciationR, extinctionR,sampleF;
extern int 			speciestreePr;
extern double 			propTemp;
extern int				curGeneration; /*prior for joint gene tree*/
extern int				activeParams[NUM_LINKED][MAX_NUM_DIVS]; /* a table holding the parameter status          */
extern int				autoClose;                              /* autoclose                                     */
extern int 				autoOverwrite;                          /* Overwrite or append outputfiles when nowarnings=yes */
extern Chain			chainParams;                            /* holds parameters for Markov chain             */
extern CharInformation	*charInfo;								/* holds critical information about characters   */
extern char				*charSetNames;                          /* holds names of character sets                 */
extern Comptree			comptreeParams;                         /* holds parameters for comparetree command      */
extern char				*constraintNames;					    /* holds names of constraints (| sep. strings)   */
extern Calibration		constraintCalibration[MAX_NUM_CONSTRAINTS];    /* holds calibration structs of constraints      */
extern int				dataType;                               /* type of data                                  */
extern int				defChars;                               /* flag for whether number of characters is known*/
extern int				defConstraints;                         /* flag for whether constraints on tree are read */
extern int				defMatrix;                              /* flag for whether matrix is successfull read   */
extern int				defPairs;                               /* flag for whether constraints on tree are read */
extern int				defPartition;                           /* flag for whether character partition is read  */
extern int				defTaxa;                                /* flag for whether number of taxa is known      */
extern Doublet			doublet[16];                            /* holds information on states for doublets      */
extern int				echoMB;							   	    /* flag used by Manual to prevent echoing        */
extern long unsigned int expecting;								/* variable denoting expected token type         */
extern int				fileNameChanged;					    /* has file name been changed?                   */
extern int				foundNewLine;                           /* whether a new line has been found             */
extern char				gapId;                                  /* gap character Id                              */
extern long int			globalSeed;                             /* seed that is initialized at start up          */
extern char				*headerNames;                           /* string to hold headers in sump and plot       */
extern int 				inComment;                              /* flag for whether input stream is commented    */
extern int				inferAncStates;					   	    /* should ancestral states be inferred (y/n)     */
extern int				inferSiteRates;					   	    /* should site rates be inferred (y/n)           */
extern int				inferPosSel;					   	    /* should positive selection be inferred (y/n)   */
extern char				inputFileName[100];                     /* input (NEXUS) file name                       */
extern int				inSumtBlock;                            /* are we in the sumt block                      */
extern int				inValidCommand;                         /* a useful flag set whenever you enter a cmd    */
extern int  			isInAmbig, isInPoly;                    /* flags whether we are within () or {}          */
extern int  			isMixed;			                    /* flags whether dataset is mixed                */
extern int				inMrbayesBlock;                         /* flag for whether we are in a mrbayes block    */
extern int				isTranslateDef;							/* is a translation block defined                */
extern int				isUserTreeDefined;                      /* flag indicating whether user tree is found    */
extern FILE				*logFileFp;                             /* file pointer to log file                      */
extern char				logFileName[100];                       /* name of the log file                          */
extern int				logToFile;                              /* should screen output be logged to a file      */
extern char				manFileName[100];						/* name of man file								 */
extern char				matchId;                                /* mach character Id                             */
extern int				*matrix;                                /* matrix containing original data               */
extern int  			matrixHasPoly;                          /* flag for whether matrix has polymorphisms     */
extern int				memAllocs[NUM_ALLOCS];                  /* allocated memory flags                        */
extern int				mode;				                    /* mode of program (interactive/noninteractive)  */
extern MoveType		    moveTypes[NUM_MOVE_TYPES];              /* holds information on the move types          */
extern char				missingId;                              /* missing character Id                          */
extern Model			modelParams[MAX_NUM_DIVS];              /* holds model params for up to 30 partitions    */
extern int				nBitsInALong;                			/* number of bits in a long			             */
extern int				noWarn;                					/* no warnings on overwriting files              */
extern int				numActiveLocks;					   	    /* number of active, locked nodes                */
extern int				numChar;                                /* number of characters in character matrix      */
extern int				numCharSets;                            /* holds number of character sets                */
extern int				numColumns;                             /* number of parameter columns for sump and plot */
extern int				numComments;                            /* number of nested comments				     */
extern int				numCurrentDivisions;                    /* number of partitions of data                  */
extern int				numDefinedConstraints;                  /* number of constraints defined                 */
extern int				numDefinedPartitions;                   /* number of partitions defined                  */
extern int				numMoveTypes;		                    /* the number of move types                     */
extern int				numOpenExeFiles;					    /* number of execute files open                  */
extern int				numPartitions;                          /* number of current partitions                  */
extern int				numRows;                                /* number of parameter rows for sump and plot    */
extern int				numTaxa;                                /* number of taxa in character matrix            */
extern int				numTaxaSets;                            /* holds number of taxa sets                     */
extern int				outGroupNum;                            /* number of outgroup taxon                      */
extern ParmInfo			paramTable[];						    /* information on parameters                     */
extern char				*partitionNames;                        /* hold names of partitions (first is "default") */
extern MrBFlt			*parameterValues;                       /* vector holding sump or plot parameters        */
extern int				partitionNum;                           /* number of current partition                   */
extern Plot				plotParams;                             /* holds parameters for plot command             */
extern int				printAncStates[MAX_NUM_DIVS];           /* divisions to print anc states for             */
extern int				quitOnError;							/* quit on error?					             */
extern int				readWord;							    /* should we read a word next?                   */
extern ReassembleInfo	reassembleParams;		                /* holds parameters for reassemble command       */
extern MrBFlt			relConstraintProbs[30];                 /* relative probs. of trees with constraint      */
extern int				replaceLogFile;                         /* should logfile be replace/appended to         */
extern long int			runIDSeed;                              /* seed used only for generating run ID [stamp]  */
extern char				spacer[10];                             /* holds blanks for printing indentations        */
extern long int			swapSeed;                               /* seed used only for determining which to swap  */
extern Sump				sumpParams;                             /* holds parameters for sump command             */
extern Sumt				sumtParams;                             /* holds parameters for sumt command             */
extern char				stamp[11];                   			/* holds a unique identifier for each analysis   */
extern char				*taxaNames;                             /* holds name of taxa                            */
extern TaxaInformation	*taxaInfo;								/* holds critical information about taxa         */
extern int				*tempSet;                               /* temporarily holds defined character set       */
extern char				*taxaSetNames;                          /* holds names of taxa sets                      */
extern int  			theAmbigChar;                           /* int containing ambiguous character            */
extern char				*transFrom;                             /* translation block information                 */
extern char				*transTo;                               /* translation block information                 */
extern int				longIntegerSize;                        /* size of an long integer                       */
extern int				userBrlensDef;                          /* are the branch lengths on user tree defined   */
extern int              userLevel;                              /* the level of the user                         */ 	
extern Tree				*userTree;                              /* user tree                                     */
extern int              userLevel;                              /* user level                                    */

#						if defined (MPI_ENABLED)
extern int 				proc_id;                                /* process ID (0, 1, ..., num_procs-1)           */
extern int 				num_procs;                              /* number of active processors                   */
extern MrBFlt			myStateInfo[4];                         /* likelihood/prior/heat vals of me              */
extern MrBFlt			partnerStateInfo[4];                    /* likelihood/prior/heat vals of partner         */
#						endif
