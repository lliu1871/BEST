

extern int      AddToString (char *s1, char *s2, int *x);
extern unsigned long int Expecting (int y);
extern int      CheckString (char *s1, char *s2, int *x);
extern int      CheckStringValidity (char *s);
extern int      DerootUserTree (TreeNode *p);
extern void     FinishTree (TreeNode *p, int *i, int isThisTreeRooted);
extern int      FreeMatrix (void);
extern int      GetNameFromString (char *s, char *tkn, int n);
extern void     GetUserDownPass (TreeNode *p, TreeNode **x, int *y);
void     		GetToken (int *tokenType);
int      		FindValidCommand (char *tk, int *numMatches);
extern int      IsArgValid (char *s, char *validArg);
int             IsIn (char ch, char *s);
int             IsSame (char *s1, char *s2);
int             IsWhite (char c);
FILE    		*OpenBinaryFileR (char *name);
FILE 			*OpenTextFileA (char *name);
FILE    		*OpenTextFileR (char *name);
FILE 			*OpenTextFileW (char *name);
int				LineTermType (FILE *fp);
int     		LongestLine (FILE *fp);
extern int      ParseCommand (char *s);
extern int      RootUserTree (TreeNode *p);
extern void     SetUpParms (void);
void     		ShowNodes (TreeNode *p, int indent, int isThisTreeRooted);
int      		ShowTree (TreeNode *r, int isThisTreeRooted, int n);
char     		WhichNuc (int x);
