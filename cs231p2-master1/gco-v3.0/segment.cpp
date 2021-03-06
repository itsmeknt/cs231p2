#include <matrix.h>
#include <mex.h>   

#include <vector>
#include "graph.h"
#include "energy.h"

using namespace std;

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  mexPrintf("1\n");

  /*
  	typedef Graph<int,int,int> GraphType;
	GraphType *g = new GraphType(2,  1); 

	g -> add_node(); 
	g -> add_node(); 

	g -> add_tweights( 0, 1, 5 );
	g -> add_tweights( 1, 2, 6 );
	g -> add_edge( 0, 1, 3, 4 );

	int flow = g -> maxflow();

	printf("Flow = %d\n", flow);
	printf("Minimum cut:\n");
	if (g->what_segment(0) == GraphType::SOURCE)
		printf("node0 is in the SOURCE set\n");
	else
		printf("node0 is in the SINK set\n");
	if (g->what_segment(1) == GraphType::SOURCE)
		printf("node1 is in the SOURCE set\n");
	else
		printf("node1 is in the SINK set\n");

	delete g;
    */

        typedef Energy<int,int,int> EnergyType;
        
  /*
  
  		// Minimize the following function of 3 binary variables:
		// E(x, y, z) = x - 2*y + 3*(1-z) - 4*x*y + 5*|y-z|
		   
		EnergyType::Var varx, vary, varz;
		EnergyType *e = new EnergyType(1,1,NULL);

		varx = e -> add_variable();
		vary = e -> add_variable();
		varz = e -> add_variable();

		e -> add_term1(varx, 0, 1);  // add term x 
		e -> add_term1(vary, 0, -2); // add term -2*y
		e -> add_term1(varz, 3, 0);  // add term 3*(1-z)

		e -> add_term2(varx, vary, 0, 0, 0, -4); // add term -4*x*y
		e -> add_term2(vary, varz, 0, 5, 5, 0); // add term 5*|y-z|

		EnergyType::TotalValue Emin = e -> minimize();
		
		printf("Minimum = %d\n", Emin);
		printf("Optimal solution:\n");
		printf("x = %d\n", e->get_var(varx));
		printf("y = %d\n", e->get_var(vary));
		printf("z = %d\n", e->get_var(varz));

		delete e;
        */
  
  int label [1];
  int N = 1;
  

  for (int i = 0; i < 1; i++) {
    // initialization
    EnergyType::Var a;
    std::vector<EnergyType::Var> vars(N);

    EnergyType *e = new EnergyType(1,1,NULL);

    // add a node
    vars[i] = e->add_variable();

    // add the unary term for a node
    e->add_term1(vars[i], 1, 1);

    // add the pairwise term for an edge
    e->add_term2(vars[i], vars[i], 1, 1, 1, 1);

    // perform energy minimization
    EnergyType::TotalValue mnE = e->minimize();

    // get new labels
    if (e->get_var(vars[i]))
      label[i] = 1;
    else
      label[i] = 0;
  }
  
  


//declare variables
    mxArray *a_in_m, *b_in_m, *c_out_m, *d_out_m;
    const mwSize *dims;
    double *a, *b, *c, *d;
    int dimx, dimy, numdims;
    int i,j;

//associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);

//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];

//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    d_out_m = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    c = mxGetPr(c_out_m);
    d = mxGetPr(d_out_m);

//do something
    for(i=0;i<dimx;i++)
    {
        for(j=0;j<dimy;j++)
        {
            mexPrintf("element[%d][%d] = %f\n",j,i,a[i*dimy+j]);
            c[i*dimy+j] = a[i*dimy+j]+5; //adds 5 to every element in a
            d[i*dimy+j] = b[i*dimy+j]*b[i*dimy+j]; //squares b
        }
    }

    return;
}
