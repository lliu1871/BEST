#include "speciestree.h"

FILE *gfopen(char *filename, char *mode)
{
   FILE *fp=(FILE*)fopen(filename, mode);
   if(fp==NULL) {
      printf("\nerror when opening file %s\n", filename);
      exit(-1);
   }
   return(fp);
}


void error2 (char * message)
{ printf("\nError: %s.\n", message); exit(-1); }

static time_t time_start;
void starttime (void)
{  
   time_start=time(NULL);
}


void strcase (char *str, int direction)
{
/* direction = 0: to lower; 1: to upper */
   char *p=str;
   if(direction)  while(*p) { *p=(char)toupper(*p); p++; }
   else           while(*p) { *p=(char)tolower(*p); p++; }
}

int matIout (FILE *fout, int x[], int n, int m)
{
   int i,j;
   for (i=0,FPN(fout); i<n; i++,FPN(fout)) 
      FOR(j,m) fprintf(fout,"%6d", x[i*m+j]);
   return (0);
}

/*///////////////////////////////////////////*/
/*      math function                       */
/*/////////////////////////////////////////*/
long factorial(int n)
{
   long f, i;
   if (n>10) error2("n>10 in factorial");
   for (i=2,f=1; i<=(long)n; i++) f*=i;
   return (f);
}


static unsigned int z_rndu=137;
static unsigned int w_rndu=123456757;

void SetSeed (unsigned int seed)
{
   if(sizeof(int)-4) puts("oh-oh, we are in trouble.  int not 32-bit?");
   z_rndu=170*(seed%178)+137;
   w_rndu = seed*127773;
}


double rndu (void)
{
/* U(0,1): AS 183: Appl. Stat. 31:188-190 
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
*/
   static unsigned int x_rndu=11, y_rndu=23;
   double r;
   

  
   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
/*
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
*/
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;

  return (r-(int)r);
}

double rndgamma1 (double s);
double rndgamma2 (double s);

double rndgamma (double s)
{
/* random standard gamma (Mean=Var=s,  with shape parameter=s, scale para=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1
*/
   double r=0;

   if (s<=0)      puts ("jgl gamma..");
   else if (s<1)  r=rndgamma1 (s);
   else if (s>1)  r=rndgamma2 (s);
   else           r=-log(rndu());
   return (r);
}


double rndgamma1 (double s)
{
/* random standard gamma for s<1
   switching method
*/
   double r, x=0,small=1e-37,w;
   static double a,p,uf,ss=10,d;

   if (s!=ss) {
      a=1-s;
      p=a/(a+s*exp(-a));
      uf=p*pow(small/a,s);
      d=a*log(a);
      ss=s;
   }
   for (;;) {
      r=rndu();
      if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
      else if (r>uf)  x=a*pow(r/p,1/s), w=x;
      else            return (0);
      r=rndu ();
      if (1-r<=w && r>0)
         if (r*(w+1)>=1 || -log(r)<=w)  continue;
      break;
   }
   return (x);
}

double rndgamma2 (double s)
{
/* random standard gamma for s>1
   Best's (1978) t distribution method
*/
   double r,d,f,g,x;
   static double b,h,ss=0;
   if (s!=ss) {
      b=s-1;
      h=sqrt(3*s-0.75);
      ss=s;
   }
   for (;;) {
      r=rndu ();
      g=r-r*r;
      f=(r-0.5)*h/sqrt(g);
      x=b+f;
      if (x <= 0) continue;
      r=rndu();
      d=64*r*r*g*g*g;
      if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
   }
   return (x);
}

int comparedouble (const void *a, const void *b)
{  
   double aa = *(double*)a, bb= *(double*)b;
   return (aa==bb ? 0 : aa > bb ? 1 : -1);
}



double Lngamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.

   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x-1;

   if((double)nx==x && nx>0 && nx<10)
      lng=log((double)factorial(nx));
   else {
      if(x<=0) {
         error2("lngamma not implemented for x<0");
         if((int)x-x==0) { puts("lngamma undefined"); return(-1); }
         for (fneg=1; x<0; x++) fneg/=x;
         if(fneg<0) error2("strange!! check lngamma");
         fneg=log(fneg);
      }
      if (x<7) {
         f=1;  z=x-1;
         while (++z<7)  f*=z;
         x=z;   f=-log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673 
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}



