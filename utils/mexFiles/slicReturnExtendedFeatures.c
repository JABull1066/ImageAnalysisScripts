//=================================================================================
//  slic.mex.c
//
//
//  AUTORIGHTS
//  Copyright (C) 2015 Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
//
//  Created by Radhakrishna Achanta on 12/01/15.
//=================================================================================
/*Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met

 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of EPFL nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

void rgbtolab(int* rin, int* gin, int* bin, int sz, double* lvec, double* avec, double* bvec)
{
    int i; int sR, sG, sB;
    double R,G,B;
    double X,Y,Z;
    double r, g, b;
    const double epsilon = 0.008856;	//actual CIE standard
    const double kappa   = 903.3;		//actual CIE standard

    const double Xr = 0.950456;	//reference white
    const double Yr = 1.0;		//reference white
    const double Zr = 1.088754;	//reference white
    double xr,yr,zr;
    double fx, fy, fz;
    double lval,aval,bval;

    for(int i = 0; i < sz; i++)
    {
        sR = rin[i]; sG = gin[i]; sB = bin[i];

        R = sR/255.0;
        G = sG/255.0;
        B = sB/255.0;

        if(R <= 0.04045)	r = R/12.92;
        else				r = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g = G/12.92;
        else				g = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b = B/12.92;
        else				b = pow((B+0.055)/1.055,2.4);

        X = r*0.4124564 + g*0.3575761 + b*0.1804375;
        Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
        Z = r*0.0193339 + g*0.1191920 + b*0.9503041;

        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;

        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;

        lval = 116.0*fy-16.0;
        aval = 500.0*(fx-fy);
        bval = 200.0*(fy-fz);

        lvec[i] = lval; avec[i] = aval; bvec[i] = bval;
    }
}

void getLABXYSeeds(int STEP, int width, int height, int* seedIndices, int* numseeds)
{
    const bool hexgrid = false;
	int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip,yerrperstrip;
    int xoff,yoff;
    int x,y;
    int xe,ye;
    int seedx,seedy;
    int i;

	xstrips = (0.5+(double)(width)/(double)(STEP));
	ystrips = (0.5+(double)(height)/(double)(STEP));

    xerr = width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = width - STEP*xstrips;}
    yerr = height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = height- STEP*ystrips;}

	xerrperstrip = (double)(xerr)/(double)(xstrips);
	yerrperstrip = (double)(yerr)/(double)(ystrips);

	xoff = STEP/2;
	yoff = STEP/2;

    n = 0;
	for( y = 0; y < ystrips; y++ )
	{
		ye = y*yerrperstrip;
		for( x = 0; x < xstrips; x++ )
		{
			xe = x*xerrperstrip;
            seedx = (x*STEP+xoff+xe);
            if(hexgrid){ seedx = x*STEP+(xoff<<(y&0x1))+xe; if(seedx >= width)seedx = width-1; }//for hex grid sampling
            seedy = (y*STEP+yoff+ye);
            i = seedy*width + seedx;
			seedIndices[n] = i;
			n++;
		}
	}
    *numseeds = n;
}

void PerformSuperpixelSLIC(double* lvec, double* avec, double* bvec, double* kseedsl, double* kseedsa, double* kseedsb, double* kseedsx, double* kseedsy, int width, int height, int numseeds, int* klabels, int STEP, double compactness)
{
    int x1, y1, x2, y2;
	double l, a, b;
	double dist;
	double distxy;
    int itr;
    int n;
    int x,y;
    int i;
    int ind;
    int r,c;
    int k;
    int sz = width*height;
	const int numk = numseeds;
	int offset = STEP;

    double* clustersize = mxMalloc(sizeof(double)*numk);
    double* inv         = mxMalloc(sizeof(double)*numk);
    double* sigmal      = mxMalloc(sizeof(double)*numk);
    double* sigmaa      = mxMalloc(sizeof(double)*numk);
    double* sigmab      = mxMalloc(sizeof(double)*numk);
    double* sigmax      = mxMalloc(sizeof(double)*numk);
    double* sigmay      = mxMalloc(sizeof(double)*numk);
    double* distvec     = mxMalloc(sizeof(double)*sz);
	double invwt = 1.0/((STEP/compactness)*(STEP/compactness));

	for( itr = 0; itr < 10; itr++ )
	{
		for(i = 0; i < sz; i++){distvec[i] = DBL_MAX;}

		for( n = 0; n < numk; n++ )
		{
            x1 = kseedsx[n]-offset; if(x1 < 0) x1 = 0;
            y1 = kseedsy[n]-offset; if(y1 < 0) y1 = 0;
            x2 = kseedsx[n]+offset; if(x2 > width)  x2 = width;
            y2 = kseedsy[n]+offset; if(y2 > height) y2 = height;

			for( y = y1; y < y2; y++ )
			{
				for( x = x1; x < x2; x++ )
				{
					i = y*width + x;

					l = lvec[i];
					a = avec[i];
					b = bvec[i];

					dist =			(l - kseedsl[n])*(l - kseedsl[n]) +
                                    (a - kseedsa[n])*(a - kseedsa[n]) +
                                    (b - kseedsb[n])*(b - kseedsb[n]);

					distxy =		(x - kseedsx[n])*(x - kseedsx[n]) + (y - kseedsy[n])*(y - kseedsy[n]);

					dist += distxy*invwt;

					if(dist < distvec[i])
					{
						distvec[i] = dist;
						klabels[i]  = n;
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
        for(k = 0; k < numk; k++)
        {
            sigmal[k] = 0;
            sigmaa[k] = 0;
            sigmab[k] = 0;
            sigmax[k] = 0;
            sigmay[k] = 0;
            clustersize[k] = 0;
        }

		ind = 0;
        for( r = 0; r < height; r++ )
        {
            for( c = 0; c < width; c++ )
            {
                if(klabels[ind] >0)
                {
                    sigmal[klabels[ind]] += lvec[ind];
                    sigmaa[klabels[ind]] += avec[ind];
                    sigmab[klabels[ind]] += bvec[ind];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    clustersize[klabels[ind]] += 1.0;
                }
                ind++;
            }
        }

		{for( k = 0; k < numk; k++ )
		{
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
		}}

		{for( k = 0; k < numk; k++ )
		{
			kseedsl[k] = sigmal[k]*inv[k];
			kseedsa[k] = sigmaa[k]*inv[k];
			kseedsb[k] = sigmab[k]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
		}}
	}
    mxFree(sigmal);
    mxFree(sigmaa);
    mxFree(sigmab);
    mxFree(sigmax);
    mxFree(sigmay);
    mxFree(clustersize);
    mxFree(inv);
    mxFree(distvec);
}

void EnforceConnectivity(int* labels, int width, int height, int numSuperpixels,int* nlabels, int* finalNumberOfLabels)
{
    int i,j,k;
    int n,c,count;
    int x,y;
    int ind;
    int label;
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    const int sz = width*height;
    int* xvec = mxMalloc(sizeof(int)*sz);
	int* yvec = mxMalloc(sizeof(int)*sz);
	const int SUPSZ = sz/numSuperpixels;
	for( i = 0; i < sz; i++ ) nlabels[i] = -1;
	int oindex = 0;
	int adjlabel = 0;//adjacent label
    label = 0;
	for( j = 0; j < height; j++ )
	{
		for( k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;

				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				{for( n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}}

				count = 1;
				for( c = 0; c < count; c++ )
				{
					for( n = 0; n < 4; n++ )
					{
						x = xvec[c] + dx4[n];
						y = yvec[c] + dy4[n];

						if( (x >= 0 && x < width) && (y >= 0 && y < height) )
						{
							int nindex = y*width + x;

							if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}

					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= SUPSZ >> 2)
				{
					for( c = 0; c < count; c++ )
					{
                        ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}

			oindex++;
		}
	}
	*finalNumberOfLabels = label;

	mxFree(xvec);
	mxFree(yvec);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 3) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if(nlhs!=7) {
        mexErrMsgIdAndTxt("SLIC:nlhs","Two outputs required, a labels and the number of labels, i.e superpixels.");
    }
    //---------------------------
    // Variable declarations
    //---------------------------
    int numSuperpixels = 200;//default value
    double compactness = 10;//default value
    int width;
    int height;
    int sz;
    int i, ii;
    int x, y, n;
    int* rin; int* gin; int* bin;
    int* klabels;
    int* clabels;
    double* lvec; double* avec; double* bvec;
    int step;
    int* seedIndices;
    int numseeds;
    double* kseedsx;double* kseedsy;
    double* kseedsl;double* kseedsa;double* kseedsb;
    double* LABMeanintensities; int* PixelCounts;
    double* SumXVector; double* SumSqrXVector; double* KVector;
    int k, kk, kkk;
    const mwSize* dims;//int* dims;
    int* outputNumSuperpixels;
    int* outlabels;
    double* outLABMeanintensities; double* outLABVariances;
    int* outPixelCounts;
    int* superpixelMinX; int* superpixelMaxX; int* superpixelMinY; int* superpixelMaxY;
    int* superpixelPerimeter;
    int* outseedsXY;
    int* neighbouringLabelsStore;
//    int* outneighbouringLabelsStore;
    int* neighbouringLabelsStoreCounter;
//    int* outneighbouringLabelsStoreCounter;
    int* twoStepneighbouringLabelsStore;
//    int* outtwoStepneighbouringLabelsStore;
    int* twoStepneighbouringLabelsStoreCounter;
//    int* outtwoStepneighbouringLabelsStoreCounter;
    double* collectedFeatures;
    double* outCollectedFeatures;
    int finalNumberOfLabels;
    unsigned char* imgbytes;
    //---------------------------
    int numelements   = mxGetNumberOfElements(prhs[0]) ;
    int numdims = mxGetNumberOfDimensions(prhs[0]) ;
    dims  = (int*)mxGetDimensions(prhs[0]) ;
    imgbytes  = (unsigned char*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    width = dims[1]; height = dims[0];//Note: first dimension provided is height and second is width
    sz = width*height;
    //---------------------------
    numSuperpixels  = mxGetScalar(prhs[1]);
    compactness     = mxGetScalar(prhs[2]);

    //---------------------------
    // Allocate memory
    //---------------------------
    rin    = mxMalloc( sizeof(int)      * sz ) ;
    gin    = mxMalloc( sizeof(int)      * sz ) ;
    bin    = mxMalloc( sizeof(int)      * sz ) ;
    //lvec    = mxMalloc( sizeof(double)      * sz ) ;
    //avec    = mxMalloc( sizeof(double)      * sz ) ;
    //bvec    = mxMalloc( sizeof(double)      * sz ) ;
    lvec    = mxCalloc( sz, sizeof(double) ) ;
    avec    = mxCalloc( sz, sizeof(double) ) ;
    bvec    = mxCalloc( sz, sizeof(double) ) ;
    klabels = mxMalloc( sizeof(int)         * sz );//original k-means labels
    clabels = mxMalloc( sizeof(int)         * sz );//corrected labels after enforcing connectivity
    seedIndices = mxMalloc( sizeof(int)     * sz );
    LABMeanintensities = mxCalloc( sz * 3,sizeof(double));
    PixelCounts = mxCalloc( sz ,sizeof(int) );
    SumXVector = mxCalloc( sz * 3,sizeof(double) ); // Used to store sum[ x_i - K ]for calculating variance
    SumSqrXVector = mxCalloc( sz * 3,sizeof(double) ); // Used to store sum of squares sum[ (x_i - K)^2 ] for calculating variance
    KVector = mxMalloc( sizeof(double)     * sz * 3); // Used to store constant K for calculating variance with shift in location to avoid catastrophic cancellation

    superpixelMinX = mxCalloc( sz ,sizeof(int));
    superpixelMinY = mxCalloc( sz ,sizeof(int));
    superpixelMaxX = mxCalloc( sz ,sizeof(int));
    superpixelMaxY = mxCalloc( sz ,sizeof(int));

    superpixelPerimeter = mxCalloc( sz ,sizeof(int));


    //---------------------------
    // Perform color conversion
    //---------------------------
    //if(2 == numdims)
    if(numelements/sz == 1)//if it is a grayscale image, copy the values directly into the lab vectors
    {
        for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                lvec[i] = imgbytes[ii];
                avec[i] = imgbytes[ii];
                bvec[i] = imgbytes[ii];
                ii++;
            }
        }
    }
    else//else covert from rgb to lab
    {
        for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                rin[i] = imgbytes[ii];
                gin[i] = imgbytes[ii+sz];
                bin[i] = imgbytes[ii+sz+sz];
                ii++;
                PixelCounts[i]=0;
            }
        }
        rgbtolab(rin,gin,bin,sz,lvec,avec,bvec);
    }


	for(x = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                if(PixelCounts[i]!=0)
                {
                	mexPrintf("Panic!");
                }
            }
        }

    //---------------------------
    // Find seeds
    //---------------------------
    step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5;
    getLABXYSeeds(step,width,height,seedIndices,&numseeds);

    kseedsx    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsy    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsl    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsa    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsb    = mxMalloc( sizeof(double)      * numseeds ) ;
    for(k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        kseedsl[k] = lvec[seedIndices[k]];
        kseedsa[k] = avec[seedIndices[k]];
        kseedsb[k] = bvec[seedIndices[k]];
    }

    //---------------------------
    // Compute superpixels
    //---------------------------
    PerformSuperpixelSLIC(lvec, avec, bvec, kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,width,height,numseeds,klabels,step,compactness);
    //---------------------------
    // Enforce connectivity
    //---------------------------
    EnforceConnectivity(klabels,width,height,numSuperpixels,clabels,&finalNumberOfLabels);
    //---------------------------
    // Assign output labels
    //---------------------------
    plhs[0] = mxCreateNumericMatrix(height,width,mxINT32_CLASS,mxREAL);
    outlabels = mxGetData(plhs[0]);

    // preparation for storing superpixel neighbours array, up to 30 per superpixel, and two step neighbours, up to 80 per superpixel
    int neighbouringLabelsStoreWidth = 30;
    neighbouringLabelsStore = mxCalloc( finalNumberOfLabels*neighbouringLabelsStoreWidth, sizeof(int) ) ;
    neighbouringLabelsStoreCounter = mxCalloc( finalNumberOfLabels, sizeof(int) ) ;
    int twoStepneighbouringLabelsStoreWidth = 80;
    twoStepneighbouringLabelsStore = mxCalloc( finalNumberOfLabels*twoStepneighbouringLabelsStoreWidth, sizeof(int) ) ;
    twoStepneighbouringLabelsStoreCounter = mxCalloc( finalNumberOfLabels, sizeof(int) ) ;


    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    for(x = 0, ii = 0; x < width; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(y = 0; y < height; y++)
        {
            i = y*width+x;
            outlabels[ii] = clabels[i];

            //------------------------------------------------------------
			// As this is the last time we loop over all pixels, once the
			// final out label has been assigned we bin the intensities
			// and use this to calculate the superpixel colour.
			//------------------------------------------------------------
			// j,k run over width and height
			// index = width*k + j

    		LABMeanintensities[3*outlabels[ii]] += lvec[i];
    		LABMeanintensities[3*outlabels[ii]+1] += avec[i];
    		LABMeanintensities[3*outlabels[ii]+2] += bvec[i];

    		// Calculate variances for LAB
    		// If entry is the first for pixel, set K to be the value for first entry
    		if(PixelCounts[clabels[i]]+1 == 1)
    		{
    			KVector[3*outlabels[ii]] = lvec[i];
    			KVector[3*outlabels[ii]+1] = avec[i];
    			KVector[3*outlabels[ii]+2] = bvec[i];

    			// Also set initial min and max X and Y for superpixel
    			superpixelMinX[clabels[i]] = x;
    			superpixelMinY[clabels[i]] = y;
    			superpixelMaxX[clabels[i]] = x;
    			superpixelMaxY[clabels[i]] = y;

    			superpixelPerimeter[clabels[i]] = 0;
    		}

    		SumXVector[3*outlabels[ii]] += (lvec[i] - KVector[3*outlabels[ii]]);// Used to store sum[ x_i - K ]for calculating variance
    		SumXVector[3*outlabels[ii]+1] += (avec[i] - KVector[3*outlabels[ii]+1]);
    		SumXVector[3*outlabels[ii]+2] += (bvec[i] - KVector[3*outlabels[ii]+2]);
    		SumSqrXVector[3*outlabels[ii]] += ((lvec[i] - KVector[3*outlabels[ii]])*(lvec[i] - KVector[3*outlabels[ii]]));
    		SumSqrXVector[3*outlabels[ii]+1] += ((avec[i] - KVector[3*outlabels[ii]+1])*(avec[i] - KVector[3*outlabels[ii]+1]));
    		SumSqrXVector[3*outlabels[ii]+2] += ((bvec[i] - KVector[3*outlabels[ii]+2])*(bvec[i] - KVector[3*outlabels[ii]+2]));

    		// While we're here, we also populate the neighbouring labels list
    		int mylabel = clabels[i];
    		int nlsStartIndex = mylabel*neighbouringLabelsStoreWidth;
    		// Find labels of surrounding pixels
    		for( n = 0; n < 4; n++ )
            {
                // Check neighbouring pixel
                int x_neighbour = x + dx4[n];
                int y_neighbour = y + dy4[n];
                // if not off the edge
                if( (x_neighbour >= 0 && x_neighbour < width) && (y_neighbour >= 0 && y_neighbour < height) )
                {
                    int neighbourindex = y_neighbour*width + x_neighbour;
                    if( clabels[neighbourindex] != mylabel )
                    {
                        // We've found a neighbour. Check if this index is in the neighbouring labels list
                        bool neighbourLabelInList = false;
                        int numNeighboursSoFar = neighbouringLabelsStoreCounter[mylabel];

                        //for( int nlsRowIndex = 0; nlsRowIndex < numNeighboursSoFar; nlsRowIndex++)
                        for( int nlsRowIndex = 0; nlsRowIndex < neighbouringLabelsStoreWidth; nlsRowIndex++)
                        {
                            if(clabels[neighbourindex] == neighbouringLabelsStore[nlsStartIndex+nlsRowIndex])
                            {
                                neighbourLabelInList = true;
                                break;
                            }
                        }
                        if(!neighbourLabelInList)
                        {
                            // We've found a neighbour we don't currently have for this superpixel! Add it to the list
                            // and increment the neighbouring label counter
                            // Also, check that we're not overflowing the list
                            assert(numNeighboursSoFar < neighbouringLabelsStoreWidth);
                            neighbouringLabelsStore[nlsStartIndex+numNeighboursSoFar] = clabels[neighbourindex];

                            neighbouringLabelsStoreCounter[mylabel]++;
                        }
                    }
                }
            }

            // Check superpixel ranges
            if (superpixelMinX[mylabel] > x)
            {
                superpixelMinX[mylabel] = x;
            }
            if (superpixelMinY[mylabel] > y)
            {
                superpixelMinY[mylabel] = y;
            }
            if (superpixelMaxX[mylabel] < x)
            {
                superpixelMaxX[mylabel] = x;
            }
            if (superpixelMaxY[mylabel] < y)
            {
                superpixelMaxY[mylabel] = y;
            }

            // Add to superpixel perimeter
            // Check neighbouring pixels: add (4 - neighbours in superpixel) to perimeter count
            int exposedEdges = 4;
            int j;
            if (x > 0)
            {
                j = y*width+(x-1);
                if (clabels[j] == clabels[i])
                {
                    exposedEdges -= 1;
                }
            }
            if (x < width)
            {
                j = y*width+(x+1);
                if (clabels[j] == clabels[i])
                {
                    exposedEdges -= 1;
                }
            }
            if (y > 0)
            {
                j = (y-1)*width+x;
                if (clabels[j] == clabels[i])
                {
                    exposedEdges -= 1;
                }
            }
            if (y < (height-1))
            {
                j = (y+1)*width+x;
                if (clabels[j] == clabels[i])
                {
                    exposedEdges -= 1;
                }
            }

            superpixelPerimeter[clabels[i]] += exposedEdges;
    		PixelCounts[clabels[i]] += 1;

            ii++;

        }
    }

    //---------------------------
    // Assign number of labels/seeds
    //---------------------------
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    outputNumSuperpixels = (int*)mxGetData(plhs[1]);//gives a void*, cast it to int*
    *outputNumSuperpixels = finalNumberOfLabels;

    //---------------------------
    // Match each superpixel with its average colour intensity for output
    //---------------------------
    plhs[2] = mxCreateNumericMatrix(3,finalNumberOfLabels,mxDOUBLE_CLASS,mxREAL);
	outLABMeanintensities = mxGetData(plhs[2]);

	plhs[3] = mxCreateNumericMatrix(1,finalNumberOfLabels,mxINT32_CLASS,mxREAL);
	outPixelCounts = mxGetData(plhs[3]);

    plhs[4] = mxCreateNumericMatrix(2,finalNumberOfLabels,mxINT32_CLASS,mxREAL);
	outseedsXY = mxGetData(plhs[4]);

    plhs[5] = mxCreateNumericMatrix(3,finalNumberOfLabels,mxDOUBLE_CLASS,mxREAL);
	outLABVariances = mxGetData(plhs[5]);

	plhs[6] = mxCreateNumericMatrix(26,finalNumberOfLabels,mxDOUBLE_CLASS,mxREAL);
	outCollectedFeatures = mxGetData(plhs[6]);




	for(k = 0; k < finalNumberOfLabels; k++)
    {
		//outseedsXY[k] = seedIndices[k];
		kk = k*2;
		kkk = k*3;
		outseedsXY[kk] = kseedsx[k];
        outseedsXY[kk+1] = kseedsy[k];
        outLABMeanintensities[kkk] = LABMeanintensities[kkk]/PixelCounts[k];
        outLABMeanintensities[kkk+1] = LABMeanintensities[kkk+1]/PixelCounts[k];
        outLABMeanintensities[kkk+2] = LABMeanintensities[kkk+2]/PixelCounts[k];
        outPixelCounts[k] = PixelCounts[k];
        outLABVariances[kkk] = (SumSqrXVector[kkk] - (SumXVector[kkk]*SumXVector[kkk]/PixelCounts[k]))/PixelCounts[k];
        outLABVariances[kkk+1] = (SumSqrXVector[kkk+1] - (SumXVector[kkk+1]*SumXVector[kkk+1]/PixelCounts[k]))/PixelCounts[k];
        outLABVariances[kkk+2] = (SumSqrXVector[kkk+2] - (SumXVector[kkk+2]*SumXVector[kkk+2]/PixelCounts[k]))/PixelCounts[k];

//        outneighbouringLabelsStoreCounter[k] = neighbouringLabelsStoreCounter[k];
        int nlsStartIndex = k*neighbouringLabelsStoreWidth;
        int nlsTwoStepStartIndex = k*twoStepneighbouringLabelsStoreWidth;
        // Loop over label's neighbours
        for(int counter = 0; counter < neighbouringLabelsStoreCounter[k]; counter++)
        {
            // Add the neighbour to the output store
            int neighbourLabel = neighbouringLabelsStore[nlsStartIndex + counter];
//            outneighbouringLabelsStore[nlsStartIndex + counter] = neighbourLabel;

            // Find the neighbouring labels of this label and loop over them!
            int neighboursNlsStartIndex = neighbourLabel*neighbouringLabelsStoreWidth;
            for(int counter2 = 0; counter2 < neighbouringLabelsStoreCounter[neighbourLabel]; counter2++)
            {
                int candidateTwoStepLabel = neighbouringLabelsStore[neighboursNlsStartIndex + counter2];
                bool labelExists = false;
                // Loop over the current twoStep store for our label, and compare with the candidate
                for(int counter3 = 0; counter3 < twoStepneighbouringLabelsStoreCounter[k]; counter3++)
                {
                    // Compare candidate with this element of list
                    if(twoStepneighbouringLabelsStore[nlsTwoStepStartIndex + counter3] == candidateTwoStepLabel)
                    {
                        labelExists = true;
                        break;
                    }
                }
                // Now, if our candidate is NOT in the list, we add it on and increase the counter
                if(!labelExists)
                {
                    twoStepneighbouringLabelsStore[nlsTwoStepStartIndex + twoStepneighbouringLabelsStoreCounter[k]] = candidateTwoStepLabel;
                    twoStepneighbouringLabelsStoreCounter[k]++;
                }

            }
        }
	}


    // OK, now that we've sorted out our one- and two- step neighbours, now we make our feature vectors!
    // 1-3 = our mean (Lab)
    // 4-6 = our variance (Lab)
    // 7-9 = one-step mean
    // 10-12 = one-step variance
    // 13-15 = two-step mean
    // 16-18 = two-step variance
    // 19-20 = superpixel width and height
    // 21 = superpixel aspect ratio (width/height)
    // 22 = number of pixels within superpixel
    collectedFeatures = mxCalloc( 26*finalNumberOfLabels, sizeof(double) ) ;
    for(k = 0; k < finalNumberOfLabels; k++)
    {
        int featuresStartIndex = 26*k;
        int meansAndVariancesStartIndex = 3*k;
        collectedFeatures[featuresStartIndex] = outLABMeanintensities[meansAndVariancesStartIndex];
        collectedFeatures[featuresStartIndex+1] = outLABMeanintensities[meansAndVariancesStartIndex+1];
        collectedFeatures[featuresStartIndex+2] = outLABMeanintensities[meansAndVariancesStartIndex+2];
        // Variances
        collectedFeatures[featuresStartIndex+3] = outLABVariances[meansAndVariancesStartIndex];
        collectedFeatures[featuresStartIndex+4] = outLABVariances[meansAndVariancesStartIndex+1];
        collectedFeatures[featuresStartIndex+5] = outLABVariances[meansAndVariancesStartIndex+2];

        int nlsStartIndex = k*neighbouringLabelsStoreWidth;
        for(int colourChannel = 0; colourChannel < 3; colourChannel++)
        {
            // Loop over one step neighbours and add their mean and variances
            // Loop over label's neighbours
            double neighbourMeansSum = 0;
            double neighbourVariancesSum = 0;
            for(int counter = 0; counter < neighbouringLabelsStoreCounter[k]; counter++)
            {
                // Add the neighbour to the output store
                int neighbourLabel = neighbouringLabelsStore[nlsStartIndex + counter];
                int neighbourMeansAndVariancesStartIndex = 3*neighbourLabel;
                neighbourMeansSum += outLABMeanintensities[neighbourMeansAndVariancesStartIndex+colourChannel];
                neighbourVariancesSum += outLABVariances[neighbourMeansAndVariancesStartIndex+colourChannel];
            }
            double neighbourMeansAverage = neighbourMeansSum/neighbouringLabelsStoreCounter[k];
            double neighbourVariancesAverage = neighbourVariancesSum/neighbouringLabelsStoreCounter[k];
            collectedFeatures[featuresStartIndex+6+colourChannel] = neighbourMeansAverage;
            collectedFeatures[featuresStartIndex+9+colourChannel] = neighbourVariancesAverage;
        }

        // Repeat for two-step neighbours
        int nlsStartIndex_twoStep = k*twoStepneighbouringLabelsStoreWidth;
        for(int colourChannel = 0; colourChannel < 3; colourChannel++)
        {
            // Loop over one step neighbours and add their mean and variances
            // Loop over label's neighbours
            double neighbourMeansSum = 0;
            double neighbourVariancesSum = 0;
            for(int counter = 0; counter < twoStepneighbouringLabelsStoreCounter[k]; counter++)
            {
                // Add the neighbour to the output store
                int neighbourLabel = twoStepneighbouringLabelsStore[nlsStartIndex_twoStep + counter];
                int neighbourMeansAndVariancesStartIndex = 3*neighbourLabel;
                neighbourMeansSum += outLABMeanintensities[neighbourMeansAndVariancesStartIndex+colourChannel];
                neighbourVariancesSum += outLABVariances[neighbourMeansAndVariancesStartIndex+colourChannel];
            }
            double neighbourMeansAverage = neighbourMeansSum/twoStepneighbouringLabelsStoreCounter[k];
            double neighbourVariancesAverage = neighbourVariancesSum/twoStepneighbouringLabelsStoreCounter[k];
            collectedFeatures[featuresStartIndex+12+colourChannel] = neighbourMeansAverage;
            collectedFeatures[featuresStartIndex+15+colourChannel] = neighbourVariancesAverage;
        }

        // Now for the superpixel size metrics
        // 19-20 = superpixel width and height
        // 21 = superpixel aspect ratio (width/height)
        // 22 = number of pixels within superpixel
        // 23 = superpixelPerimeter
        // 24 = perimeter / area ratio
        collectedFeatures[featuresStartIndex+18] = superpixelMaxX[k] - superpixelMinX[k];
        collectedFeatures[featuresStartIndex+19] = superpixelMaxY[k] - superpixelMinY[k];
        collectedFeatures[featuresStartIndex+20] = 1.0*(superpixelMaxX[k] - superpixelMinX[k])/(superpixelMaxY[k] - superpixelMinY[k]);
        collectedFeatures[featuresStartIndex+21] = PixelCounts[k];
        collectedFeatures[featuresStartIndex+22] = superpixelPerimeter[k];
        collectedFeatures[featuresStartIndex+23] = 1.0*superpixelPerimeter[k] / PixelCounts[k];

        // 25-26 = number of neighbours / two-step neighbours
        collectedFeatures[featuresStartIndex+24] = neighbouringLabelsStoreCounter[k];
        collectedFeatures[featuresStartIndex+25] = twoStepneighbouringLabelsStoreCounter[k];


        // Now that we've sorted the collected features for this label, we update the output list
        for(int featureCounter = 0; featureCounter < 26; featureCounter++)
        {
            outCollectedFeatures[featuresStartIndex+featureCounter] = collectedFeatures[featuresStartIndex+featureCounter];
        }
    }


    //---------------------------
    // Deallocate memory
    //---------------------------
    mxFree(rin);
    mxFree(gin);
    mxFree(bin);
    mxFree(lvec);
    mxFree(avec);
    mxFree(bvec);
    mxFree(klabels);
    mxFree(clabels);
    mxFree(seedIndices);
    mxFree(PixelCounts);
    mxFree(LABMeanintensities);
    mxFree(kseedsx);
    mxFree(kseedsy);
    mxFree(kseedsl);
    mxFree(kseedsa);
    mxFree(kseedsb);
    mxFree(SumXVector);
    mxFree(SumSqrXVector);
    mxFree(KVector);
    mxFree(neighbouringLabelsStore);
    mxFree(neighbouringLabelsStoreCounter);
    mxFree(twoStepneighbouringLabelsStore);
    mxFree(twoStepneighbouringLabelsStoreCounter);
    mxFree(collectedFeatures);
    mxFree(superpixelMinX);
    mxFree(superpixelMinY);
    mxFree(superpixelMaxX);
    mxFree(superpixelMaxY);
    mxFree(superpixelPerimeter);
}

//PixelCounts
