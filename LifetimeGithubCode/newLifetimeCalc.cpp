//NEW lifetime calculation to include spectral smoothing
// PROGRAM THAT CALCULATES THE LIFETIME OF ELECTRONS
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>
using namespace std;

//typedef std::pair<double, double> Point;//defining a "point" for the Ramer-Douglas-Peucker algorithm

int numGraphs = 1000; //number of graphs to be averaged
int n;
double cut = 300.;
double tauelec_preampB = 91.5031*1e-6; // PRE_AMP FOR CATHODE (ASSUMING AMPS ARE "SWITCHED")
double tauelec_preampA = 43.4835*1e-6; //PRE_AMP FOR ANODE (ASSUMING AMPS ARE "SWITCHED")
double gain_AoverB = 0.89512; 
int nBeginCh4 = 4400; //Beginning of cathode minimum
int nBeginCh3 = 7800; //Beginning of anode maximum
int sampleSize = 3000;
double epsilon = 0.00001;
int varHalfWidth = 20;
int noiseSampleSize = 4000;
double numStDev = 1.0;
int beginFit1 = 100;
int endFit1 = 300;
int beginFit2 = 300;
int endFit2 = 500;
int beginFit3 = 500;
int endFit3 = 730;
int beginFit4 = 730;
int endFit4 = 900;
int fitBufferSize = 30;
int numPointsForMeanOfEndPoint = 10;


double getCorrection(double fittedDriftTime, double tauelec, double taulife){
  double laura = (1.-exp(-(fittedDriftTime)*(1./tauelec + 1./taulife)))/(fittedDriftTime*(1./tauelec + 1./taulife));
  return laura;
  
}

double variance(double *array, double mean, int size){
  double var = 0;
  for (int i=0; i<size; i++){
    var += (array[i]-mean)*(array[i]-mean)/(size-1);
  } 
  return var;
}

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

double t2Method3(double *yValues, double *xValues, int nMax){
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
  double a = (sumy*sumx2 - sumx*sumxy)/denom;
  double b = (npoints*sumxy - sumx*sumy)/denom;
  return -a/b;
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

double lifetime(double Qa, double Qc, double t1, double t2, double t3){
  double tau = (t2+(t1+t3)/2)/log(fabs(Qa/Qc));
  cout<<"Qa = "<<Qa<<" Qc = "<<Qc<<" t1 = "<<t1<<" t2 = "<<t2<<" t3 = "<<t3<<" lifetime = "<<tau<<endl;
  return tau;
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
  fitGraphReal->Fit(reFit4, "CQMR+");  
  fitGraphImaginary->Fit(imFit1, "CQR");
  fitGraphImaginary->Fit(imFit2, "CQR+");
  fitGraphImaginary->Fit(imFit3, "CQR+");
  fitGraphImaginary->Fit(imFit4, "CQMR+");
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
        if (j == 0){
        }
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
  double stdOfVar = sqrt(variance(varArray, meanOfVar, 1000));
  cout<<stdOfVar<<endl;
  double tolerance = meanOfVar+(2.0*stdOfVar);
  cout<<tolerance<<endl;
  int count = 0;
  int startingN;
  for (int i=0; i<arraySize; i++){
    enhancedYValues[i] = yValues[i];
    if(varArray[i]>(tolerance)){
      if (count == 0){
        startingN = i;
      }
      count++;
    }
  }
  double startVal = yValues[startingN-1];
  double endVal;
  for (int i=0; i<numPointsForMeanOfEndPoint; i++){
    endVal += yValues[startingN+count+i]/numPointsForMeanOfEndPoint;
  }
  for (int i=0; i<count; i++){
    enhancedYValues[i+startingN] = yValues[startingN+i-(count+1)] + (endVal-startVal)*i/count;
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


void newLifetimeCalc(){
  TFile *ch3 = new TFile("/home/gmattimoe/data/K-890GK-800GA800A1000_50.100.200V.cm_FibreIn_23.08.ch3.traces.root", "read");
  TFile *ch4 = new TFile("/home/gmattimoe/data/K-890GK-800GA800A1000_50.100.200V.cm_FibreIn_23.08.ch4.traces.root", "read");

  TGraph *g1 = (TGraph*)ch3->Get("graph1"); //importing graph
  n = g1->GetN(); //number of points in graph
  double *x = g1->GetX(); // time values

  double *totalCh3 = new double[n]();
  totalCh3 = totalGraphs(ch3, totalCh3);
  double medianCh3 = findMedian(totalCh3, sampleSize);
  for (int i=0; i<n; i++){
    totalCh3[i] = totalCh3[i] - medianCh3;
  }

  double *totalCh4 = new double[n]();
  totalCh4 = totalGraphs(ch4, totalCh4);
  double medianCh4 = findMedian(totalCh4, sampleSize);
  for (int i=0; i<n; i++){
    totalCh4[i] = totalCh4[i] - medianCh4;
  }  
  
  double *varianceArrayCh3 = new double[n]();
  localVariance(totalCh3, n, varianceArrayCh3, varHalfWidth);
  TGraph *varianceGraphCh3 = new TGraph(n, x, varianceArrayCh3);  
  double *enhancedTotalCh3 = new double[n]();
  zeroFluctuationRemovalCh3(totalCh3, enhancedTotalCh3, varianceArrayCh3, n); 
  
  double *varianceArrayCh4 = new double[n]();
  localVariance(totalCh4, n, varianceArrayCh4, varHalfWidth);
  TGraph *varianceGraphCh4 = new TGraph(n, x, varianceArrayCh4);
  double *enhancedTotalCh4 = new double[n]();
  zeroFluctuationRemoval(totalCh4, enhancedTotalCh4, varianceArrayCh4, n);
  
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
  
  int maxnCh4 = maxNValue(fourierRemovalCh4, nBeginCh4, 500);
  int maxnCh3 = maxNValue(fourierRemovalCh3, nBeginCh3, 500);
  double Qc = fourierRemovalCh4[maxnCh4];
  double Qa = fourierRemovalCh3[maxnCh3];
  cout<<"Qc = "<<Qc<<endl;
  Qa /= gain_AoverB; //ASSUMING THE AMPS ARE "SWITCHED"
  int count = 0;
  double newQa;
  double newQc;
  double correctionA;
  double correctionC;
  double t1 = x[maxnCh4];
  double t2 = t2Method3(fourierRemovalCh3, x, maxnCh3) - t1;
  double t3 = x[maxnCh3];
  double tau = lifetime(Qa, Qc, t1, t2, t3); //initial estimate of lifetime
  while (count<20){
    correctionA = getCorrection(t3, tauelec_preampB, tau);
    correctionC = getCorrection(t1, tauelec_preampA, tau);
    newQa = Qa/correctionA;
    newQc = Qc/correctionC;
    tau = lifetime(newQa, newQc, t1, t2, t3);      
    count++; 
  }
  cout<<"lifetime = "<<tau<<endl;
  
  
  TCanvas *c2 = new TCanvas("c1", "multipads", 900,700);
  TGraph *anodeGraph = new TGraph(n, x, fourierRemovalCh3);
  TGraph *cathodeGraph = new TGraph(n, x, fourierRemovalCh4);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(anodeGraph);
  mg1->Add(cathodeGraph);
  mg1->Draw("ALP");
}