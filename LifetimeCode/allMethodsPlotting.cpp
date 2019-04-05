#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>
using namespace std;

//typedef std::pair<double, double> Point;//defining a "point" for the Ramer-Douglas-Peucker algorithm

int numGraphs = 1000; //number of graphs to be averaged
int aveBound = 10; //number of points either side of the point in the moving average
int n;
int nBeginCh4 = 4400; //Beginning of cathode minimum
int nBeginCh3 = 7800; //Beginning of anode maximum
int nSubPoints = 500; //number of points to search for average over
int nAve = 3; //number of values to be averaged either side of nMax
int sampleSize = 3000;
double epsilon = 0.00001;
double cut = 200.;
int varHalfWidth = 20;
int noiseSampleSize = 4000;
double numStDev = 1.0;
int beginFit1 = 100;
int endFit1 = 300;
int beginFit2 = 300;
int endFit2 = 500;
int beginFit3 = 500;
int endFit3 = 720;
int beginFit4 = 720;
int endFit4 = 900;
int fitBufferSize = 10;
int numPointsForMeanOfEndPoint = 10;



int maxNValue(double *yValues, int nBegin, int nSubPoints){
  double max = 0;
  int maxn;
  for (int i=0; i<nSubPoints; i++){
    if (std::fabs(yValues[i+nBegin]) > std::abs(max)){
      max = yValues[i+nBegin];
      maxn = i+nBegin;
    }
  }
  return maxn;
}

double t2Method1(double *yValues, double *xValues, int nMax){
  int ipoint = nMax;
  double y1 = yValues[ipoint];
  double y = yValues[ipoint-1]; 
  while (y < y1){
    ipoint--;
    y1 = y;
    y = yValues[ipoint-1];
  }
  return xValues[ipoint];
}

double t2Method2(double *yValues, double *xValues, int nMax){
  double max = yValues[nMax];
  int ipoint = nMax;
  while (yValues[ipoint] > 0.1*max){
    ipoint--;
  }
  double x0 = xValues[ipoint]+(xValues[ipoint]-xValues[nMax])*yValues[ipoint]/(yValues[nMax]-yValues[ipoint]); 
  return x0;
}

double t2Method3(double *yValues, double *xValues, int nMax, double *a, double *b){
  double max = yValues[nMax];
  int ipoint = nMax;
  while (yValues[ipoint] > 0.1*max){
    ipoint--;
  }
  double xi;
  double yi;
  double sumx = 0.0;
  double sumy = 0.0;
  double sumxy = 0.0;
  double sumx2 = 0.0;
  int npoints = nMax+1-ipoint;
  for (int i=0; i<npoints; i++){
    xi = xValues[ipoint+i];
    yi = yValues[ipoint+i];
    sumx += xi;
    sumy += yi;
    sumxy += xi*yi;
    sumx2 += xi*xi;
  }
  double denom = (npoints*sumx2 - sumx*sumx);
  *a = (sumy*sumx2 - sumx*sumxy)/denom;
  *b = (npoints*sumxy - sumx*sumy)/denom;
  return -(*a)/(*b);
}


double maxValueAve(int nMax, double *yValues){
  double maxAve = 0;
  for (int i=-nAve; i<=nAve; i++){
    maxAve += yValues[i+nMax]/(2*nAve+1);
  }
  return maxAve;
}


double *movingAverage(double *original, double *movingAve){
  for (int i=0;i<n;i++){
    if (i<aveBound || i>=n-aveBound){
      movingAve[i]=original[i];
    }
    else{
      movingAve[i]=0;
      for (int j=-aveBound;j<=aveBound; j++){
	      movingAve[i] += original[i+j]/(2*aveBound+1);
      }
    }
  }
  return movingAve;
}


double *totalGraphs(TFile *file, double *total){
  for (int igraph=1; igraph<=numGraphs; igraph++){
    string name = Form("graph%i", igraph);
    TGraph *g = (TGraph*)file->Get(name.c_str());
    double *y = g->GetY();
    for (int ipoint = 0; ipoint<n; ipoint++){
      total[ipoint] += y[ipoint];
    }
  }
  for (int i=0; i<n; i++){
    total[i] = total[i]/numGraphs;
  }
  return total;
}


double logistic(int x, double y){ //function to 
  if (x <= cut){
    return y;
  }
  else{
    return 0.0;
  }  
}


double mean(double *array, int size){
  double sum = 0;
  for (int i=0; i<size; i++){
    sum += array[i];
  } 
  double mean = sum/size;
  return mean;
}


double findMedian(double *array, int sizeOfArray) { 
    const int size = sizeOfArray;
    int index[size];
    TMath::Sort(size, array, index, false); 
    if (size % 2 != 0){
       return array[index[size-1/2]]; 
    }  
    else{
      return (array[index[(size-2)/2]] + array[index[size/2]])/2.0; 
    }
}


double variance(double *array, double mean, int size){
  double var = 0;
  for (int i=0; i<size; i++){
    var += (array[i]-mean)*(array[i]-mean)/(size-1);
  } 
  return var;
}


double *fourierConv(double *original, double *fourier){
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &n, "r2c");
  fftr2c->SetPoints(original);
  fftr2c->Transform();  
  double *re = new double[n];
  double *im = new double[n];
  fftr2c->GetPointsComplex(re,im);
  for (int i=0; i<n; i++){
    re[i] = logistic(i,re[i]);
    im[i] = logistic(i,im[i]);
  }
  TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &n, "c2r");  
  fftc2r->SetPointsComplex(re,im);
  fftc2r->Transform();
  fftc2r->GetPoints(fourier);
  for (int i=0; i<n; i++){
    fourier[i] = fourier[i]/n;
  }
  return fourier;
}
 
 
double PerpendicularDistance(const pair<double, double> &pt, const pair<double, double> &lineStart, const pair<double, double> &lineEnd){
  double nx = (-lineEnd.second + lineStart.second);
  double ny = (lineEnd.first - lineStart.first);
  double dist = fabs(((pt.first - lineStart.first)*nx + (pt.second - lineStart.second)*ny)/pow((nx*nx + ny*ny),0.5));
  return dist;
}

 
void RamerDouglasPeucker(vector<pair<double,double>> &pointList, double epsilon, vector<pair<double,double>> &out){
	if(pointList.size()<2)
		throw invalid_argument("Not enough points to simplify");
 
	// Find the point with the maximum distance from line between start and end
	double dmax = 0.0;
	size_t index = 0;
	size_t end = pointList.size()-1;
	for(size_t i = 1; i < end; i++)
	{
		double d = PerpendicularDistance(pointList[i], pointList[0], pointList[end]);
		if (d > dmax)
		{
			index = i;
			dmax = d;
		}
	}
 
	// If max distance is greater than epsilon, recursively simplify
	if(dmax > epsilon)
	{
		// Recursive call
		vector<pair<double, double>> recResults1;
		vector<pair<double, double>> recResults2;
		vector<pair<double, double>> firstLine(pointList.begin(), pointList.begin()+index+1);
		vector<pair<double, double>> lastLine(pointList.begin()+index, pointList.end());
		RamerDouglasPeucker(firstLine, epsilon, recResults1);
		RamerDouglasPeucker(lastLine, epsilon, recResults2);
 
		// Build the result list
		out.assign(recResults1.begin(), recResults1.end()-1);
		out.insert(out.end(), recResults2.begin(), recResults2.end());
		if(out.size()<2)
			throw runtime_error("Problem assembling output");
	} 
	else 
	{
		//Just return start and end points
		out.clear();
		out.push_back(pointList[0]);
		out.push_back(pointList[end]);
	}
}

void justFourierTransform(double *original, double *re, double *im, int inputSize){
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &inputSize, "r2c");
  fftr2c->SetPoints(original);
  fftr2c->Transform();
  fftr2c->GetPointsComplex(re,im);
}

void inverseFourierTransform(double *re, double *im, double *newTimeDomainArray, int outputSize){
  int size = outputSize;
  TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &size, "c2r");  
  fftc2r->SetPointsComplex(re,im);
  fftc2r->Transform();
  fftc2r->GetPoints(newTimeDomainArray);
  for (int i=0; i<size; i++){
    newTimeDomainArray[i] = newTimeDomainArray[i]/size;
  }
}

void spectralRemoval(double *signalSpectrum, double *noiseSpectrum, TGraph *removedSpectralSignal, int signalArraySize, int noiseArraySize, double noiseTolerance){
  int noiseInt;
  int noiseCount;
  int pointInt;
  int signalCount;
  double newSpectrumValue;
  double startValue;
  double endValue;
  for (int i=0; i<signalArraySize; i++){
    pointInt = removedSpectralSignal->GetN();
    noiseCount = 0;
    signalCount = 0;
    noiseInt = i*noiseArraySize/signalArraySize;
    if (i<endFit4 && i>beginFit1){
      while (noiseSpectrum[noiseInt] > noiseTolerance){
        noiseCount++;
        noiseInt++;
      }
      if (noiseCount>0){
        signalCount = noiseCount*signalArraySize/noiseArraySize;
        i += (signalCount-1);
      }
      else{
        removedSpectralSignal->SetPoint(pointInt, i, signalSpectrum[i]);
      }
    }  
  }
}

void fourierRemovalConvolution(double *originalSignal, double *noiseSignal, double *finalSmoothedSignal, int signalArraySize, int noiseArraySize){
  int specSize = n/2 + 1;
  int noiseSpecSize = noiseSampleSize/2 + 1;
  double *re = new double[specSize]();
  double *im = new double[specSize]();
  double *reNoise = new double[noiseSpecSize]();
  double *imNoise = new double[noiseSpecSize]();
  justFourierTransform(originalSignal, re, im, signalArraySize);
  justFourierTransform(noiseSignal, reNoise, imNoise, noiseArraySize);
  double noiseMeanRe = 0.0;
  double noiseMeanIm = 0.0;
  for (int i=0; i<specSize; i++){
    noiseMeanRe += reNoise[i]/specSize;
    noiseMeanIm += imNoise[i]/specSize;
  }
  double reNoiseSTD = sqrt(variance(reNoise, noiseMeanRe, noiseSpecSize));
  double imNoiseSTD = sqrt(variance(imNoise, noiseMeanIm, noiseSpecSize));
  double noiseToleranceRe = noiseMeanRe + numStDev*reNoiseSTD;
  double noiseToleranceIm = noiseMeanIm + numStDev*imNoiseSTD;  
  TGraph *fitGraphReal = new TGraph();
  TGraph *fitGraphImaginary = new TGraph();
  spectralRemoval(re, reNoise, fitGraphReal, specSize, noiseSpecSize, reNoiseSTD);
  spectralRemoval(im, imNoise, fitGraphImaginary, specSize, noiseSpecSize, imNoiseSTD);
  TF1 *reFit1 = new TF1("reFit1", "pol6", beginFit1, endFit1+fitBufferSize);
  TF1 *reFit2 = new TF1("reFit2", "pol6", beginFit2, endFit2+fitBufferSize);
  TF1 *reFit3 = new TF1("reFit3", "pol6", beginFit3, endFit3+fitBufferSize);
  TF1 *reFit4 = new TF1("reFit4", "pol6", beginFit4, endFit4+fitBufferSize);
  TF1 *imFit1 = new TF1("imFit1", "pol6", beginFit1, endFit1+fitBufferSize);
  TF1 *imFit2 = new TF1("imFit2", "pol6", beginFit2, endFit2+fitBufferSize);
  TF1 *imFit3 = new TF1("imFit3", "pol6", beginFit3, endFit3+fitBufferSize);
  TF1 *imFit4 = new TF1("imFit4", "pol6", beginFit4, endFit4+fitBufferSize);  
  fitGraphReal->Fit(reFit1, "CQR");
  fitGraphReal->Fit(reFit2, "CQR+");
  fitGraphReal->Fit(reFit3, "CQR+");
  fitGraphReal->Fit(reFit4, "CQR+");  
  fitGraphImaginary->Fit(imFit1, "CQR");
  fitGraphImaginary->Fit(imFit2, "CQR+");
  fitGraphImaginary->Fit(imFit3, "CQR+");
  fitGraphImaginary->Fit(imFit4, "CQR+");
  double *reNew = new double[specSize]();
  double *imNew = new double[specSize]();
  for (int i=0; i<specSize; i++){
    if (i<=beginFit1){
      reNew[i] = re[i];
      imNew[i] = im[i];
    }
    else if(i<=endFit1){
      reNew[i] = reFit1->Eval(i);
      imNew[i] = imFit1->Eval(i);
    }
    else if(i<=endFit2){
      reNew[i] = reFit2->Eval(i);
      imNew[i] = imFit2->Eval(i);
    }
    else if(i<=endFit3){
      reNew[i] = reFit3->Eval(i);
      imNew[i] = imFit3->Eval(i);
    }
    else if(i<=endFit4){
      reNew[i] = reFit4->Eval(i);
      imNew[i] = imFit4->Eval(i);
    }
    else{
      reNew[i] = 0.0;
      imNew[i] = 0.0;
    }
  }
  inverseFourierTransform(reNew, imNew, finalSmoothedSignal, signalArraySize);
}

void localVariance(double *yValues, int ySize, double *varArray, int halfWidthOfVariance){
  double *localArray = new double[2*halfWidthOfVariance+1]();
  double mean;
  double var;
  for(int i=0; i<ySize; i++){
    if (i<halfWidthOfVariance || i>=(ySize-halfWidthOfVariance)){
      varArray[i] = 0.0;
    }
    else{
      mean = 0.0;
      for(int j=-halfWidthOfVariance; j<=halfWidthOfVariance; j++){
        localArray[j+halfWidthOfVariance] = yValues[i+j];
        mean += yValues[i+j]/(2*halfWidthOfVariance+1);
      }
      var = variance(localArray, mean, (2*halfWidthOfVariance+1));
      varArray[i] = var;
    }
  }
}

void zeroFluctuationRemoval(double *yValues, double *enhancedYValues, double *varArray, int arraySize){
  double meanOfVar = 0.0;
  for (int i=0; i<arraySize; i++){
    meanOfVar += varArray[i]/n;
  }
  double stdOfVar = sqrt(variance(varArray, meanOfVar, n));
  cout<<stdOfVar<<endl;
  double tolerance = meanOfVar+(3.0*stdOfVar);
  cout<<tolerance<<endl;
  int count = 0;
  int startingN;
  for (int i=0; i<arraySize; i++){
    enhancedYValues[i] = yValues[i];
    if(varArray[i]>(tolerance)){
      while(varArray[i]>meanOfVar){
        if (count == 0){
          startingN = i;
        }
        count++;
        i++;
      }
    }
  }
  double startVal = yValues[startingN-1];
  double endVal=yValues[startingN+count];//0.0;
  for (int i=0; i<numPointsForMeanOfEndPoint; i++){
    endVal += yValues[startingN+count+i]/numPointsForMeanOfEndPoint;
  }
  for (int i=0; i<=count; i++){
    enhancedYValues[i+startingN] = yValues[startingN+i-(count+10)] + (endVal-startVal)*i/count;
  }
}

void zeroFluctuationRemovalCh3(double *yValues, double *enhancedYValues, double *varArray, int arraySize){
  double meanOfVar = 0.0;
  for (int i=0; i<arraySize; i++){
    meanOfVar += varArray[i]/n;
  }
  double stdOfVar = sqrt(variance(varArray, meanOfVar, n));
  cout<<stdOfVar<<endl;
  double tolerance = meanOfVar+(1.0*stdOfVar);
  cout<<tolerance<<endl;
  int count = 0;
  int startingN;
  for (int i=0; i<arraySize; i++){
    enhancedYValues[i] = yValues[i];
    if(varArray[i]>(tolerance)){
      while(varArray[i]>meanOfVar){
        if (count == 0){
          startingN = i;
        }
        count++;
        i++;
      }
    }
  }
  for (int i=0; i<=count; i++){
    enhancedYValues[i+startingN] = 0.0;
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN PROGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





void allMethodsPlotting(){
  //importing the root directory 
  TFile *ch3 = new TFile("/home/gmattimoe/data/K-890GK-800GA800A1000_50.100.200V.cm_FibreIn_23.08.ch3.traces.root", "read");
  TFile *ch4 = new TFile("/home/gmattimoe/data/K-890GK-800GA800A1000_50.100.200V.cm_FibreIn_23.08.ch4.traces.root", "read");

  TGraph *g1 = (TGraph*)ch3->Get("graph1"); //importing graph
  n = g1->GetN(); //number of points in graph
  double *x = g1->GetX(); // time values

  double *totalCh3 = new double[n]();
  totalCh3 = totalGraphs(ch3, totalCh3);  
  //double meanCh3 = mean(totalCh3, sampleSize);
  //double sdCh3 = sqrt(variance(totalCh3, meanCh3, sampleSize));
  double medianCh3 = findMedian(totalCh3, sampleSize);
  
  for (int i=0; i<n; i++){
    x[i] *= 1000.;
    totalCh3[i] = totalCh3[i] - medianCh3;
  }
  cout<<"CH3"<<endl;
  double *varianceArrayCh3 = new double[n]();
  localVariance(totalCh3, n, varianceArrayCh3, varHalfWidth);
  TGraph *varianceGraphCh3 = new TGraph(n, x, varianceArrayCh3);  
  double *enhancedTotalCh3 = new double[n]();
  zeroFluctuationRemovalCh3(totalCh3, enhancedTotalCh3, varianceArrayCh3, n);  
 
  double *movingAveCh3 = new double[n]();
  movingAveCh3 = movingAverage(totalCh3, movingAveCh3);
  double *fourierCh3 = new double [n]();
  fourierCh3 = fourierConv(totalCh3, fourierCh3);
  double leastSqrA;
  double leastSqrB;
  int maxn = maxNValue(fourierCh3, nBeginCh3, nSubPoints);
  double t2one = t2Method1(fourierCh3, x, maxn);
  double t2two = t2Method2(fourierCh3, x, maxn);
  double t2three = t2Method3(fourierCh3, x, maxn, &leastSqrA, &leastSqrB);
//~~~~~~~~~~~~~Ramer�Douglas�Peucker algorithm~~~~~~~~~~~~~~~~~~~~~
	// /*
  vector<pair<double,double>> pointListCh3;
	vector<pair<double,double>> pointListOutCh3;
  for (int i=0; i<n; i++){
    pointListCh3.push_back(pair<double, double> (x[i], fourierCh3[i]));
  }
  RamerDouglasPeucker(pointListCh3, epsilon, pointListOutCh3);
  
  int rdpCh3Size = pointListOutCh3.size();
  double *rdpCh3x = new double[rdpCh3Size]();
  double *rdpCh3y = new double[rdpCh3Size]();
  for (int i=0; i<rdpCh3Size; i++){
    rdpCh3x[i] = pointListOutCh3[i].first;
    rdpCh3y[i] = pointListOutCh3[i].second;
  }
  // */

//~~~~~~~~~~~~~~~Graph Least Squares Fit~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  int leastSqrN = 1000;
  int leastSqrNBegin = 7750;
  double *leastSqrX = new double[leastSqrN]();
  double *leastSqrY = new double[leastSqrN]();
  for (int i=0; i<leastSqrN; i++){
    leastSqrX[i] = x[leastSqrNBegin+i];
    leastSqrY[i] = x[leastSqrNBegin+i]*leastSqrB + leastSqrA;
  }

  double *totalCh4 = new double[n]();
  totalCh4 = totalGraphs(ch4, totalCh4);  
  //double meanCh4 = mean(totalCh4, sampleSize);
  //double sdCh4 = sqrt(variance(totalCh4, meanCh4, sampleSize));
  double medianCh4 = findMedian(totalCh4, sampleSize);
  
  for (int i=0; i<n; i++){
    totalCh4[i] = totalCh4[i] - medianCh4;
  }
  
  double *varianceArrayCh4 = new double[n]();
  localVariance(totalCh4, n, varianceArrayCh4, varHalfWidth);
  TGraph *varianceGraphCh4 = new TGraph(n, x, varianceArrayCh4);
  double *enhancedTotalCh4 = new double[n]();
  zeroFluctuationRemoval(totalCh4, enhancedTotalCh4, varianceArrayCh4, n);

  
  double *movingAveCh4 = new double[n]();
  movingAveCh4 = movingAverage(totalCh4, movingAveCh4);
  double *fourierCh4 = new double[n]();
  fourierCh4 = fourierConv(totalCh4, fourierCh4);
  
//~~~~~~~~~~~~~Ramer�Douglas�Peucker algorithm~~~~~~~~~~~~~~~~~~~~~
  // /*
	vector<pair<double,double>> pointListCh4;
	vector<pair<double,double>> pointListOutCh4;
  for (int i=0; i<n; i++){
    pointListCh4.push_back(pair<double, double> (x[i], fourierCh4[i]));
  }
  RamerDouglasPeucker(pointListCh4, epsilon, pointListOutCh4);
  
  int rdpCh4Size = pointListOutCh4.size();
  double *rdpCh4x = new double[rdpCh4Size]();
  double *rdpCh4y = new double[rdpCh4Size]();
  for (int i=0; i<rdpCh4Size; i++){
    rdpCh4x[i] = pointListOutCh4[i].first;
    rdpCh4y[i] = pointListOutCh4[i].second;
  }
  // */


//~~~~~~~~~~~~~~~Spectral Smoothing~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double *noiseTotalCh3 = new double[noiseSampleSize]();
  double *noiseTotalCh4 = new double[noiseSampleSize]();
  for (int i=0; i<noiseSampleSize; i++){
    noiseTotalCh3[i] = totalCh3[i+(n-1)-noiseSampleSize];
    noiseTotalCh4[i] = totalCh4[i+(n-1)-noiseSampleSize];
  } 
  double *fourierRemovalCh3 = new double[n]();
  fourierRemovalConvolution(enhancedTotalCh3, noiseTotalCh3, fourierRemovalCh3, n, noiseSampleSize);
  TGraph *smoothedGraphCh3 = new TGraph(n, x, fourierRemovalCh3);
  
  double *fourierRemovalCh4 = new double[n]();
  fourierRemovalConvolution(enhancedTotalCh4, noiseTotalCh4, fourierRemovalCh4, n, noiseSampleSize);
  TGraph *smoothedGraphCh4 = new TGraph(n, x, fourierRemovalCh4);  
  
//~~~~~~~~~~~~~~~Graph Formatting Anode~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  TCanvas *c2 = new TCanvas("c1", "multipads", 1300,500);
  //c2->Divide(1,2);
   
  TGraph *totCh3 = new TGraph(n,x,totalCh3);
  TGraph *movAveCh3 = new TGraph(n,x,movingAveCh3);
  TGraph *fourierGraphCh3 = new TGraph(n,x,fourierCh3); 
  TGraph *rdpCh3 = new TGraph(rdpCh3Size, rdpCh3x, rdpCh3y);
  TGraph *leastSqr = new TGraph(leastSqrN,leastSqrX,leastSqrY);

  //totCh3->SetMarkerSize(0.01);
  totCh3->SetMarkerStyle(1);
  totCh3->SetMarkerColor(9);
  totCh3->SetLineColor(9);
  totCh3->SetTitle("Anode");
  
  
  smoothedGraphCh3->SetMarkerStyle(1);
  smoothedGraphCh3->SetMarkerColor(2);
  smoothedGraphCh3->SetLineColor(1);
  smoothedGraphCh3->SetLineWidth(2);

  smoothedGraphCh4->SetMarkerStyle(1);
  smoothedGraphCh4->SetMarkerColor(2);
  smoothedGraphCh4->SetLineColor(1);
  smoothedGraphCh4->SetLineWidth(2);
  
  
  movAveCh3->SetMarkerSize(0.01);
  movAveCh3->SetMarkerStyle(1);
  movAveCh3->SetMarkerColor(1);
  movAveCh3->SetTitle("Anode Moving Mean");

  //fourierGraphCh3->SetMarkerSize(0.01);
  fourierGraphCh3->SetMarkerStyle(1);
  fourierGraphCh3->SetMarkerColor(2);
  fourierGraphCh3->SetTitle("Anode Fourier Convolution");
  
  rdpCh3->SetMarkerSize(1);
  rdpCh3->SetDrawOption("L");
  rdpCh3->SetMarkerStyle(1);
  rdpCh3->SetMarkerColor(1);
  rdpCh3->SetTitle("Anode RDP algo");
  
  leastSqr->SetLineColor(2);
  leastSqr->SetLineStyle(7);
  leastSqr->SetLineWidth(2);
  leastSqr->SetTitle("Least Squares Fit");
    
  delete[] totalCh3;
  delete[] movingAveCh3;
  delete[] fourierCh3;
  delete[] rdpCh3x;
  delete[] rdpCh3y;

//~~~~~~~~~~~~~Graph Formatting Cathode~~~~~~~~~~~~~~~~~~~~~~~~~~

  TGraph *movAveCh4 = new TGraph(n,x,movingAveCh4);
  TGraph *totCh4 = new TGraph(n,x,totalCh4);
  TGraph *enhancedTotCh4Graph = new TGraph(n,x,enhancedTotalCh4);  
  TGraph *fourierGraphCh4 = new TGraph(n,x,fourierCh4);
  TGraph *rdpCh4 = new TGraph(rdpCh4Size, rdpCh4x, rdpCh4y);
  
  //TF1 *interferenceFit1 = new TF1("interferenceFit1", "landau", beginInterferenceFit1, endInterferenceFit1);
  //totCh4->Fit(interferenceFit1, "RM");
  //totCh4->SetMarkerSize(0.01);
  totCh4->SetMarkerStyle(1);
  totCh4->SetMarkerColor(9);
  totCh4->SetLineColor(9);
  totCh4->SetTitle("Cathode");
  
  enhancedTotCh4Graph->SetMarkerStyle(1);
  enhancedTotCh4Graph->SetMarkerColor(9);
  enhancedTotCh4Graph->SetLineColor(4);
  enhancedTotCh4Graph->SetTitle("Cathode");  
  
  movAveCh4->SetMarkerSize(0.01);
  movAveCh4->SetMarkerStyle(1);
  movAveCh4->SetMarkerColor(1);
  movAveCh4->SetTitle("Cathode Moving Mean");  
  
  //fourierGraphCh4->SetMarkerSize(0.01);
  fourierGraphCh4->SetMarkerStyle(1);
  fourierGraphCh4->SetMarkerColor(2);
  fourierGraphCh4->SetTitle("Cathode Fourier Convolution"); 
  
  rdpCh4->SetMarkerSize(0.01);
  rdpCh4->SetMarkerStyle(1);
  rdpCh4->SetMarkerColor(1);
  rdpCh4->SetTitle("Cathode RDP algo");   

//~~~~~~~~~~~~~~~~~~Graph Plotting~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //c2->cd(1);
  TMultiGraph *mg1 = new TMultiGraph();
  //mg1->SetTitle("Single Waveform of Anode Signal");
  mg1->SetTitle("Mean of Cathode Waveforms with Smoothed Spectrum");
  //mg1->SetTitle("Anode Pulse Showing Least Squares Fit Pulse");
  mg1->Add(totCh3, "L");
  mg1->Add(totCh4, "L");
  //mg1->Add(enhancedTotCh4Graph, "L");
  //mg1->Add(fourierGraphCh3, "P");
  //mg1->Add(fourierGraphCh4, "P");
  mg1->Add(smoothedGraphCh3, "L");
  mg1->Add(smoothedGraphCh4, "l");
  //mg1->Add(movAveCh3, "L");
  //mg1->Add(movAveCh4, "L");
  //mg1->Add(leastSqr, "L");
  mg1->Draw("ALP");
  mg1->GetXaxis()->SetTitle("Time ms");
  mg1->GetYaxis()->SetTitle("Voltage V");
  //mg1->GetXaxis()->SetLimits(0.3,0.55);
  mg1->GetYaxis()->SetRangeUser(-0.04, 0.05);
  gPad->SetGrid(1, 1); 
  gPad->Update();
  
  
  auto legend1 = new TLegend(0.65,0.75,0.9,0.9);
  legend1->SetTextSize(0.03);
  legend1->AddEntry(totCh3, "Anode Average Signal", "l");
  legend1->AddEntry(totCh4, "Cathode Average Signal", "l");
  //legend1->AddEntry(leastSqr, "Least Squares Linear Fit", "l");
  legend1->AddEntry(smoothedGraphCh3, "Anode Smoothed Signal", "l");
  legend1->AddEntry(smoothedGraphCh4, "Cathode Smoothed Signal", "l");
  legend1->Draw(); 
  
  //gPad->BuildLegend();
  
  /*c2->cd(2);
  TMultiGraph *mg2 = new TMultiGraph();
  //mg2->Add(varianceGraphCh4, "L");
  mg2->Add(enhancedTotCh4Graph, "L");
  //mg2->Add(smoothedGraphCh4, "P");
  //mg2->Add(varianceGraphCh3, "L");
  mg2->SetTitle("Mean of 1000 waveforms of the Cathode Signal With Interference Removed");
  //mg2->Add(totCh3);
  //mg2->Add(totCh4);
  //mg2->Add(fourierGraphCh3);
  //mg2->Add(fourierGraphCh4);
  //mg2->Add(movAveCh3);
  //mg2->Add(movAveCh4);
  mg2->Draw("APL");
  mg2->GetXaxis()->SetTitle("Time ms");
  mg2->GetYaxis()->SetTitle("Voltage V");
  mg2->GetYaxis()->SetRangeUser(-0.04, 0.01);
  //mg2->GetYaxis()->SetRangeUser(0.0, 0.00005);
  mg2->GetXaxis()->SetLimits(-0.1,0.2);
  gPad->SetGrid(1, 1); 
  gPad->Update();   
  //gPad->Modified();
  //gPad->BuildLegend();
  */
  delete [] movingAveCh4;
  delete [] totalCh4;

}