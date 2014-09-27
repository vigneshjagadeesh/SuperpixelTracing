// Computing Optical Flow through a mex wrapper to opencv ... 
// Basic Idea is to get the template and target frames, followed by flow
// computation

// Assuming there are only single channel double images ... 

// My Include Path - C:\OpenCV2.3\build\include
// Lib Path        - C:\OpenCV2.3\build\x64\vc9\lib
// mex 

#include "mex.h"

#include <iostream>
#include <fstream>

#include "opencv2/core/core.hpp"
#include "opencv2/video/tracking.hpp"

using namespace std;
using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Variable Initializations
    double *inp_param;
    unsigned char *srcImg, *dstImg;
    int i;
    if ( (nrhs != 3)  ) 
        mexErrMsgTxt("ERROR1");
    
    // Set up all input pointers 
    inp_param = (double*)mxGetPr(prhs[0]);
    srcImg = (unsigned char*)mxGetPr(prhs[1]);    
    dstImg = (unsigned char*)mxGetPr(prhs[2]);
    
    // Set up all the input parameters ....
    int numRows = (int)inp_param[0];
    int numCols = (int)inp_param[1];
    double pyrScale = inp_param[2];
    int levels      = (int)inp_param[3];
    int winSize     = (int)inp_param[4];
    int iterations  = (int)inp_param[5];
    int polyN       = (int)inp_param[6];
    double polySigma= inp_param[7];
    mexPrintf("NumRows %d \nNumCols %d\n pyrScale %f\n levels %d\n winSize %d\n iterations %d\n polyN %d\n polySigma %f\n", numRows, numCols, pyrScale, levels, winSize, iterations, polyN, polySigma);
    Mat I = Mat::zeros(numRows, numCols, CV_8UC1);
    Mat J = Mat::zeros(numRows, numCols, CV_8UC1);
    Mat flow = Mat::zeros(numRows, numCols, CV_32FC2);
    
    for( int rowIter = 0; rowIter < numRows; rowIter++)
        for( int colIter = 0; colIter < numCols; colIter++)
        {
            I.at<unsigned char> (rowIter, colIter) = srcImg[numRows*colIter + rowIter];
            J.at<unsigned char> (rowIter, colIter) = dstImg[numRows*colIter + rowIter];
        } 
    
    calcOpticalFlowFarneback(I, J, flow, pyrScale, levels, winSize, iterations, polyN, polySigma, OPTFLOW_USE_INITIAL_FLOW );
        
    // RETURN VARIABLE
    int numFlow = numRows*numCols;
    plhs[0] = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
    double* outArray1 = (double*)mxGetPr(plhs[0]);
    double* outArray2 = (double*)mxGetPr(plhs[1]);
    for( int rowIter = 0; rowIter < numRows; rowIter++)
    {
        float* flowPtr = flow.ptr<float> (rowIter);
        for( int colIter = 0; colIter < numCols; colIter++)
        {
            outArray1[numRows*colIter + rowIter] = (double) flowPtr[2*colIter];
            outArray2[numRows*colIter + rowIter] = (double) flowPtr[2*colIter + 1];            
        }
    }
    return;
}