#include "mex.h"
#include "opencv2/core/core.hpp"

#include <iostream>
#include <fstream>
using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int *inp_param;
    unsigned char *srcImg, *dstImg;
    int i;
    if ( (nrhs != 3)  ) 
        mexErrMsgTxt("ERROR1");
    
    // Set up all input pointers 
    inp_param = (int*)mxGetPr(prhs[0]);
    srcImg = (unsigned char*)mxGetPr(prhs[1]);    
    dstImg = (unsigned char*)mxGetPr(prhs[2]);
    
    // Set up all the input parameters ....
    int numRows = inp_param[0];
    int numCols = inp_param[1];
    double pyrScale = inp_param[2];
    int levels      = inp_param[3];
    int winSize     = inp_param[4];
    int iterations  = inp_param[5];
    int polyN       = inp_param[6];
    double polySigma= inp_param[7];
    
    cv::Mat I = cv::Mat::zeros(10, 10, CV_8UC1); 
    cout << "Hello World" << endl;
}