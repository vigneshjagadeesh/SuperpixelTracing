// Mex Function Interfacing with Kohli et al's robust higher order potentials code ...

#include "expand.h"
#include <stdio.h>
#include <stdlib.h>
#include "conio.h"
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    srand(40000);
        
    mxArray    *tmp;
    int        ifield, nfields,rowLen, colLen;
    int        index_iter, superpixel_iterator,i,j;
    mwIndex    jstruct;
    mwSize     NStructElems;
    mwSize     ndim;
    double     *xValues;
    double     *t_links, *n_weights;
    double     *n_links;
    int        *inp_param;
    
    // Superpixel Assignments ...
    inp_param = (int*)mxGetPr(prhs[0]);
     t_links =   (double*)mxGetPr(prhs[1]);
     n_links =   (double*)mxGetPr(prhs[2]);
     n_weights = (double*)mxGetPr(prhs[3]);
     nfields =    mxGetNumberOfFields(prhs[4]); // will be two for superpixels
     NStructElems = mxGetNumberOfElements(prhs[4]); 
     
    Energy<double> *energy = new Energy<double>(inp_param[0], inp_param[1], inp_param[2], NStructElems);
    printf("No_Objects:%d, No_Nodes:%d, No_nlinks:%d, No_SPixels:%d", inp_param[0], inp_param[1], inp_param[2], NStructElems);
    
    // Assigning Unary Potential Costs ...
    for(i = 0; i < energy->nvar*energy->nlabel; i++) 
	{
		energy->unaryCost[i] = t_links[i];
        //printf("\n %f\n", energy->unaryCost[i]);
	}
    
       
    //    //initialize grid topology - horizontal connections
    
    for(i = 0;i < inp_param[2] ;i++)
    {
        energy->pairIndex[2*i]   = (int)n_links[2*i];
        energy->pairIndex[2*i+1] = (int)n_links[2*i + 1];
        energy->pairCost[i] = n_weights[i];
        //printf("StartVertex:%d, End Vertex:%d, Weights:%f \n", (int)n_links[2*i], (int)n_links[2*i + 1], n_weights[i]);
    }

    
    
    ifield = 0;
    superpixel_iterator = 0;
        for(jstruct = 0; jstruct < NStructElems; jstruct++) {   
            tmp = mxGetFieldByNumber(prhs[4], jstruct, ifield); 
            xValues = mxGetPr(tmp);
            energy->higherElements[superpixel_iterator] = (int)xValues[0];
            //mexPrintf("\n%d \n",(int)energy->higherElements[superpixel_iterator]);            
            superpixel_iterator = superpixel_iterator + 1;
        }
    energy->AllocateHigherIndexes();
    
   
    // Assigning the indices of pixels present in every superpixel
    ifield = 1;
    superpixel_iterator = 0;
        for(jstruct = 0; jstruct < NStructElems; jstruct++) {   
            tmp = mxGetFieldByNumber(prhs[4], jstruct, ifield); 
            rowLen = mxGetN(tmp);             colLen = mxGetM(tmp);
            //mexPrintf("Row Len: %d and Col Len: %d \n",rowLen,colLen);
            
            xValues = mxGetPr(tmp);        
            for(index_iter=0; index_iter<(rowLen*colLen); index_iter++) {
                energy->higherIndex[superpixel_iterator][index_iter] = (int)xValues[index_iter];
                //mexPrintf("\n%d iter %d contains %d\n",superpixel_iterator,index_iter ,energy->higherIndex[superpixel_iterator][index_iter]);
            }            
        superpixel_iterator = superpixel_iterator + 1;
        }
    
    
    
    
// Actual Higher Order Clique Reduction ... 
//initialize truncation ratio Q, gamma_k and gamma_max for each clique
    
	for(i = 0; i < energy->nhigher; i++)
	{
        ifield = 2;
        tmp = mxGetFieldByNumber(prhs[4], i, ifield); 
        xValues = mxGetPr(tmp);
        //printf("\n tempo: %d \n", i);
		//truncation ratio 30%
		energy->higherTruncation[i] = xValues[energy->nlabel + 1];
        
		//gamma_k
		for(j = 0; j < energy->nlabel; j++) energy->higherCost[i * (energy->nlabel + 1) + j] = xValues[j]; //(rand() % 100) * 0.1;

        //gamma_max
		energy->higherCost[i * (energy->nlabel + 1) + energy->nlabel] = xValues[energy->nlabel]; //100 * 0.1;
	}

// //initialize alpha expansion - max 10 iterations
 	AExpand *expand = new AExpand(energy, 10);

//initialize solution
	int *solution = new int[inp_param[1]];
	memset(solution, 0, inp_param[1] * sizeof(int));
// 
// //solve CRF
 	expand->minimize(solution);

//print solution
    plhs[0] = mxCreateDoubleMatrix(inp_param[1], 1, mxREAL);
    double* outArray = (double*)mxGetPr(plhs[0]);
    for(i = 0; i < inp_param[1]; i++)
	{
		  outArray[i] = (double)solution[i]; 
	}
    

 
// //free memory
 	delete[] solution;
 	delete expand;
    delete energy;
    mexPrintf("Done Deal \n");
}

