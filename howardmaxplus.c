/* Howard Policy Iteration Algorithm for Max Plus Matrices

 Dec 22, 1997
 Stephane.Gaubert@inria.fr

 The function is documented in the manual page Howard.help

 The max-plus version of Howard's algorithm is described
 in the paper: 
  Numerical computation of spectral elements in max-plus algebra
            Jean Cochet-Terrasson, Guy Cohen, Stephane Gaubert,
	    Michael Mc Gettrick, Jean-Pierre Quadrat
  (to appear in the IFAC Workshop on System Structure and Control,
  Nantes, July 1997).
 */

/* 

   Howard will terminate with an error if MAX_NIterations 
   occur.
   (on all pratical random examples, even for matrices of large size
   - we tried for matrices of size <64Mb -,
   the algorithm requires "at most" 80 iterations).

*/

/* how to define infvalue and epsilon*/

#define MAX_NIterations 1000
#define INFVALUE 1000000000000000

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "mex.h"
/*#include "Display2.h"*/



#define EPSILON -1000000000000000
/*#define EPSILON  -1.79769313486231570e308*/




/* GLOBAL VARIABLES */

static int *ij;
static double *a;
static double *t;
static int nnodes;
static int narcs;
static double *chi;
static double *v;
static int *pi;
static int *NIterations;
static int *NComponents;
static int verbosemode;


/****************************************************************/

static int *newpi=NULL; /*  new policy */
/* the inverse policy is coded by a linearly chained list */
static 	int *piinv_idx=NULL; /* piinv_idx[i]= pointer to the chain of inverses
				of node i */
static	int *piinv_succ=NULL;/* piinv_idx[j]= pointer to the next inverse */
static	int *piinv_elem=NULL;/* corresponding  node  */
static	int *piinv_last=NULL;/* piinv_last[i]= last inverse of i */
static	double *c=NULL;
static  double *tau=NULL;
static  double *newtau=NULL;
static	double *vaux=NULL;
static	double *newc=NULL;
static	double *newchi=NULL;
static 	int *visited=NULL;
static	int *component=NULL;
static  double lambda=0;
static  double epsilon=0;
static  int color=1;


/* The termination tests are performed up to an epsilon
constant, which is fixed heuristically by the following routine
*/

void Epsilon(double *a,int narcs,double *epsilon)
{ 
  int i;
  double maxx,minn;
  maxx=a[0];
  minn=a[0];
  for (i=1;i<narcs; i++)
    {
      if (a[i]> maxx)
	{
	  maxx=a[i];
	}
      if (a[i]<minn)
	{
	  minn=a[i];
	}
    }
  /* this value is somehow arbitrary. */
   *epsilon=(maxx-minn+1)*0.000000001;
   if (verbosemode>1)
     {
       printf("MAX of coefs =%f MIN of coefs=%f MAX-MIN= %f epsilon=%f\n",maxx,minn,maxx-minn, *epsilon);
     }
}

/* Build an admissible policy pi and its associated cost vector c */
/* from ij  and A */

/* reasonable greedy rule to determine the first policy. 
  pi(node i) = arc with maximal weigth starting from i 
  for full random matrices, this choice of initial policy seems
  to cut the number of iterations by a factor 1.5, by comparison
  with a random initial policy */ 

static void Initial_Policy()
{ int i;  
  /* we loose a O(nnodes) time here ... */
  /* we use the auxiliary variable vaux to compute the row max of A */ 
  for (i=0; i< nnodes; i++) 
    {
      vaux[i] = EPSILON;
    }
  for (i=0; i<narcs; i++)
    {
      if (vaux[ij[i*2]] <= a[i]) 
	{
       pi[ij[i*2]] =  ij[i*2+1];
       c[ij[i*2]] = a[i];
       vaux[ij[i*2]] = a[i];
     }
    }
}



static void New_Display_Inverse()
{ 
 int i,a;
 printf("DISPLAYING INVERSE nnodes=%d\n",nnodes);
 for (i=0; i<nnodes; i++)
   { a=piinv_idx[i];
     while (a != -1)
       {
	 printf("inverse of %d= %d\n",i,piinv_elem[a]);
	 a=piinv_succ[a];
       }
   }
}

static void New_Build_Inverse()
{ int i,j,locus;
  int ptr=0;
 if (verbosemode >1 )
    {
      printf("BUILDING INVERSE\n");
    }
 for (i=0;i<nnodes;i++)
   {
     piinv_idx[i] =-1;
     piinv_last[i]=-1;
   }
 for (i=0;i<nnodes;i++)
   {
     j=pi[i];
     if (piinv_idx[j]==-1) 
       {
	 piinv_succ[ptr]=-1;
	 piinv_elem[ptr]=i;
	 piinv_last[j]=ptr;
	 piinv_idx[j]=ptr;
	 ptr++;
       }
     else
       {
       piinv_succ[ptr]=-1;
       piinv_elem[ptr]=i;
       locus=piinv_last[j];
       piinv_succ[locus]=ptr;
       piinv_last[j]=ptr;
       ptr++;
       };
   }
 if (verbosemode >1)
   {
     New_Display_Inverse();
     printf("INVERSE OK\n");
   }
}

static void Init_Depth_First()
{
  int j;
  for (j=0;j<nnodes;j++)
    {
      visited[j]=0;
      component[j]=0;
    }
}


/* Given the value of v at initial point i, we compute v[j] */
/* for all predecessor j of i, according to the spectral equation */
/* v[j]+ lambda = A(arc from j to i) + v[i] */
/* the array visited is changed by side effect */


static void New_Depth_First_Label(int i)
{ 
  int nexti,a;
  a=piinv_idx[i];
  while (a != -1 && visited[piinv_elem[a]]==0)
    {
      nexti=piinv_elem[a];
      visited[nexti]=1;
      v[nexti]= -lambda + c[nexti]+ v[i];
      component[nexti]=color;
      chi[nexti]= lambda;
      New_Depth_First_Label(nexti);
      a=piinv_succ[a];
    }
}



static void Visit_From(int initialpoint,int color)
{ int index,newindex,i;
  double weight;
  int length;
  if (verbosemode>1)
    {
      printf("visiting from node %d color=%d\n",initialpoint, color);
    }
  index=initialpoint;
  component[index]=color;
  newindex=pi[index];
  while (component[newindex]==0)
    {
      component[newindex]=color;
      index=newindex;
      newindex=pi[index];
    }
  /* a cycle has been detected, since newindex is already visited */ 
  weight=0;
  length=0;
  i=index;
  do
    {
      weight+=c[i];
      length++;
      i=pi[i];
    }
  while (i !=index);
  lambda=weight/length;
  v[i]=vaux[i]; /* keeping the previous value */
  chi[i]=lambda;
  New_Depth_First_Label(index);
}





/*       Value() */
/* Computes the value (v,chi) associated with a policy pi */

static void Value()
{ int initialpoint;
  color=1;
  if (verbosemode>0)
    {
      printf("Computing the value\n");
    }
  Init_Depth_First();
  initialpoint=0;
  do
    { 
      Visit_From(initialpoint,color);
      while ((initialpoint<nnodes) && (component[initialpoint] !=0))
	{ 
	  initialpoint++;
	}
      color++;
    }
  while (initialpoint<nnodes);
  *NComponents=--color;
  if (verbosemode>0)
    {
      printf("Value OK\n");
    }
}





static void Show_Info_Improve_Chi(int i)
{
int I,J;
/*printf("type 1 improvement\n");*/
if (verbosemode>0)
  {
  I=ij[i*2];
  J=ij[i*2+1];
  printf("Improvement of the cycle time at node %d\n",I);
  printf("arc %d: %d--->%d chi[%d]-chi[%d]=%f-%f=%g>0\n",i,I,J,I,J,chi[J],chi[I],chi[J]-chi[I]);
  }
}

static void Show_Info_Improve_Bias(int i)
{
int I,J;
/*printf("type 2 improvement\n");*/
if (verbosemode>0)
  {
  I=ij[i*2];
  J=ij[i*2+1];
  printf("Improvement of the BIAS at node %d\n",I);
  printf("A[%d]+v[%d] - chi[%d]-v[%d]= %f + %f -%f -%f =%f >0\n",i,J,I,I,a[i],v[J],chi[I],v[I],a[i]+v[J]-chi[I]-v[I]);
  }
}



static void Init_Improve()
{
  int i;
  for (i=0;i<nnodes; i++)
    {
      newchi[i]=chi[i];
      vaux[i]=v[i];
      newpi[i]=pi[i];
      newc[i]=c[i];
    }
}



static void First_Order_Improvement(int *improved)
{ int i;
  for (i=0;i<narcs; i++)
    {
      if (chi[ij[i*2+1]]>newchi[ij[i*2]])
	{ 
	  Show_Info_Improve_Chi(i);
	  *improved=1;
	  newpi[ij[i*2]]=ij[i*2+1];
	  newchi[ij[i*2]]=chi[ij[i*2+1]];
	  newc[ij[i*2]]=a[i];
	}
    }
}


static void Second_Order_Improvement(int *improved)
{ int i;     
  double w;
 if (*NComponents >1) /* a bit more complicated */
   {
     for (i=0;i<narcs; i++)
	{
	  if (chi[ij[i*2+1]]==newchi[ij[i*2]])
	  /* arc i is critical */
	{
	  w=a[i]+ v[ij[i*2+1]] - chi[ij[i*2]];
	  if (w>vaux[ij[i*2]] + epsilon)
	    {
	      Show_Info_Improve_Bias(i);
	      *improved=1;
	      vaux[ij[i*2]]=w;
	      newpi[ij[i*2]]=ij[i*2+1];
	      newc[ij[i*2]]=a[i];
	    }
	}
	}
   }
 else /* we know that all the arcs realize the max in the 
	 first order improvement */
   {
     for (i=0;i<narcs; i++)
	{
	  w=a[i]+ v[ij[i*2+1]] - chi[ij[i*2]];
	  if (w>vaux[ij[i*2]] + epsilon)
	    {
	      Show_Info_Improve_Bias(i);
	      *improved=1;
	      vaux[ij[i*2]]=w;
	      newpi[ij[i*2]]=ij[i*2+1];
	      newc[ij[i*2]]=a[i];
	    }
	}
   }
}



static void Improve(int *improved)
{
  *improved=0;
  Init_Improve();
  if (*NComponents>1) /* a first order policy improvement may occur */
    First_Order_Improvement(improved);
  if (*improved ==0)
    {
      Second_Order_Improvement(improved);
    }
}




static void Allocate_Memory()
{
  newpi=(int *)calloc(nnodes, sizeof(int));
  piinv_idx=(int *)calloc(nnodes, sizeof(int));
  piinv_succ=(int *)calloc(nnodes, sizeof(int));
  piinv_elem=(int *)calloc(nnodes, sizeof(int));
  piinv_last=(int *)calloc(nnodes, sizeof(int));
  visited=(int *)calloc(nnodes, sizeof(int));
  component=(int *)calloc(nnodes, sizeof(int));
  c=(double *)calloc(nnodes, sizeof(double));
  newc=(double *)calloc(nnodes, sizeof(double));
  vaux=(double *)calloc(nnodes, sizeof(double));
  newchi=(double *)calloc(nnodes, sizeof(double));
  if ((newchi==NULL)|| (vaux==NULL) || (newc==NULL)|| (c==NULL)|| (component==NULL) || (visited==NULL)|| (piinv_idx==NULL) ||(piinv_succ== NULL) ||( piinv_elem==NULL)|| (  piinv_last==NULL) || (newpi==NULL))
    {
      printf("error in Howard... memory allocation failed...\n");
      exit(1);
    }
}



static void Free_Memory()
{
  free(newpi);
  free(piinv_idx);
  free(piinv_succ);
  free(piinv_elem);
  free(piinv_last);
  free(visited);
  free(component);
  free(c);
  free(newc);
  free(vaux);
  free(newchi);
}










static int Check_Rows()
{ 
  int i;
  int *u=NULL;
  u=(int *)calloc(nnodes,sizeof(int));
  for (i=0; i<narcs;i++ )
    {
       u[ij[2*i]]=1;
    }
  for (i=0; i<nnodes;i++ )
    {
     if (u[i]==0)
       {
	 if (verbosemode>0)
	   {
	     printf("ERROR : node numbered %d has no sucessor (recall that nodes are internally numbered from 0)\n",i);
	   }
	 free(u);
	 return(1);
       }
    }
 free(u);
 return(0);
}

static int Universal_Security_Check()
{
   if (nnodes<1)
    {
      if (verbosemode>0)
	{
	printf("ERROR: nnodes must be a positive integer\n");
	}
      return(3);
    }
  if (narcs<1)
    { if (verbosemode>0)
      {
	printf("ERROR: narcs must be a positive integer\n");
      }
      return(4);
    }
  if (verbosemode >-1)
    {
      return(Check_Rows());
    }
  else
     {
       return(0);
     }
 }

static int Security_Check()
{ int errorflag;
  errorflag=Universal_Security_Check();
  /*if (verbosemode >0)
     {
       Display_Sparse_Matrix(narcs,ij,a);

     }*/
  return(errorflag);
}



void Import_Arguments(int *IJ, double *A,int NNODES, int NARCS, double *CHI, double *V, int *POLICY, int *NITERATIONS, int *NCOMPONENTS, int VERBOSEMODE)
{
ij=IJ;
a=A;
nnodes=NNODES;
narcs=NARCS;
chi=CHI;
v=V;
pi=POLICY;
NIterations=NITERATIONS;
NComponents=NCOMPONENTS;
verbosemode=VERBOSEMODE;
}



static void Update_Policy()
{
register int i;
      for (i=0;i<nnodes;i++)
	   {
	     pi[i]=newpi[i];
    	     c[i]=newc[i];
	     vaux[i]=v[i]; /* We need a to keep a copy of the current value function */
	   }
}



static int End_Message()
{
  if (*NIterations ==MAX_NIterations)
    {
      printf("ERROR : maximal number of iterations (=%d) reached\n",MAX_NIterations);
      printf("This should not happen (usual number of iterations is <80,\n");
      printf("the worst case number of iterations is expected to be of order nnodes\n");
      printf("If you do expect a large number of iterations, you may change\n");
      printf("manually the constant MAX_NIterations in HOWARD?.c\n");
      printf("If you suspect an error, you may wish to rerun the program with the flag verbosemode >0 to trace what happens\n");
      printf("Bug reports should be sent to Stephane.Gaubert@inria.fr\n");
      return(5);
    }
  if (verbosemode >1)
     {
       printf("END OF HOWARD: OK\n");
     }
  return(0);
}



/* INPUT of Howard Algorithm =
   ij,A,nnodes,narcs : sparse description of a matrix 
   verbosemode: normal value is zero
   higher values (1,2) yield more info for debugging purposes.
   Value -1 of verbosemode suppresses the checking of consistency 
   of data.
   
   OUTPUT =
   chi cycle time vector
   v bias
   pi optimal policy
   NIterations: Number of iterations of the algorithm
   NComponents: Number of connected components of the optimal policy

   REQUIRES: O(nnodes) SPACE
   One iteration requires: O(narcs+nnodes) TIME
   
   Experimentally, the number of iterations N_H(nnodes,...) seems to grow
   slowy with the dimension, something like 
   N_H(nnodes)=O(log(nnodes)) for full matrices.
   The matrix A must have at least one finite entry
   per row (this is checked when verbosemode is >=1)

   */


/*  
The following variables should be defined in the environment
from which the Howard routine is called.

INPUT VARIABLES
int NNODES;   
              number of nodes of the graph 
int NARCS;    
              number of arcs of the graph 
int *IJ; 
              array of integers of size 2*narcs 
       
	      for (0 <=k <narcs), the arc numbered k  goes from 
	      IJ[k][0] =(IJ[2k]) to IJ[k][1] (=IJ[2k+1])

double *A; 
              array of double of size narcs
	      A[k]=weight of the arc numbered k 

OUTPUT VARIABLES

double *V;
               array of double of size nnodes (the bias vector)
double *CHI;
               array of double of size nnodes (the cycle time vector)
int *POLICY;       
               array of integer of size nnodes (an optimal policy)
int NITERATIONS; 
               integer: the number of iterations of the algorithm

int NCOMPONENTS; 
               integer: the number of connected components of the optimal
	       policy which is returned.
*/

static int Howard(int *IJ, double *A,int NNODES,int NARCS,double *CHI,double *V,int *POLICY,int *NITERATIONS,int *NCOMPONENTS,int VERBOSEMODE)
{
  int improved=0;
  int ERRORTYPE;
  *NITERATIONS=0;
  Import_Arguments(IJ,A,NNODES,NARCS,CHI,V,POLICY,NITERATIONS,NCOMPONENTS,VERBOSEMODE);
  ERRORTYPE=Security_Check();
  if (ERRORTYPE !=0)
    {
      return(ERRORTYPE);
    }
  Allocate_Memory();
  Epsilon(a,narcs,&epsilon);
  Initial_Policy();
  New_Build_Inverse();
  do
    { 
       Value();
       /*Show_Info();*/
       Improve(&improved);
       Update_Policy();
       New_Build_Inverse();
       (*NIterations)++;
     }
  while ((improved != 0) && *NIterations <MAX_NIterations);
  Free_Memory();
  ERRORTYPE=End_Message();
  return(ERRORTYPE);
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
/* input */
double *a1;

/*output*/
double *chi1;
double *v1;
    
/*other variables used in Howard*/    

int *pi1=NULL; /* policy: dimension=nnodes */
int nnodes1;

int *ij1=NULL; /* description of arcs*/
double *amod=NULL; /* weight of arcs */
int narcs1=0, arcindex=0; /*narcs=dimension of the above*/

int niters; /* do not need extra work*/
int ncomps;
int Errortype1;
    
int i,j;    
    
    /*double *yp; 
    double *t,*y; */
    mwSize m,n; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 1) { 
	mexErrMsgTxt("One input matrix required.");
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments.");
    } 
    
    /* Check the input */ 
    
    
    
    m = mxGetM(prhs[0]); 
    n= mxGetN(prhs[0]);
    
    
    if (mxIsComplex(prhs[0]) ||  (m!= n)) 
    { mexErrMsgTxt("Program requires SQUARE matrix as input argument."); 
    }
    
    a1=mxGetPr(prhs[0]);
    
    nnodes1=n;
    
    pi1=(int *)calloc(nnodes1, sizeof(int));
    
    if (pi1==NULL)
    { 
         mexErrMsgTxt("Could not allocate memory for pi.");
    }
    
    for (i=0;i<n;i++)
        for(j=0;j<n;j++)
            if (a1[n*i+j]>-INFVALUE)
            {narcs1=narcs1+1;}
    
    ij1=(int *)calloc(2*narcs1, sizeof(int));
    amod=(double *)calloc(narcs1, sizeof(double));
    
    if ((ij1==NULL)||(amod==NULL))
    {
        mexErrMsgTxt("Could not allocate memory for ij and amod."); 
    }
    
    
    
    for (i=0;i<n;i++)
        for(j=0;j<n;j++)
            if (a1[n*i+j]>-INFVALUE)
            {amod[arcindex]=a1[n*i+j];
             ij1[2*arcindex]=j;
             ij1[2*arcindex+1]=i;
             arcindex=arcindex+1;}
    
    
    /* The double vectors of chi and v*/ 
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    
    
    
    /* Assign pointers to the various parameters */ 
    chi1=mxGetPr(plhs[0]);
    v1=mxGetPr(plhs[1]);
    
    /* Do the actual computations in a subroutine */
    Errortype1=Howard(ij1, amod, nnodes1, narcs1, chi1, v1, pi1, 
            &niters, &ncomps, 0);
    
     
    free(pi1);
    free(ij1);
    free(amod);
    return;
    
}






