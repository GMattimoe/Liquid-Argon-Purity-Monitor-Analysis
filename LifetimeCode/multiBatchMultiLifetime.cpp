// Extension of multiLifetime.C where that function is looped over all batch sizes

// THIS PROGRAM CALULATES MULTIPLE GRAPHS FROM A SAMPLE OF THE 1000 GRAPHS

// WITH THE T2 CALCULATION THE MAX VALUE WAS CALCULATED USING THE MAX VALUE FROM THE
// FOURIER CONVOLUTION DATA SET AS THE NOISE IS MUCH MORE PRONOUNCED AND THE MAX VALUE 
// IS SOMETIMES FAR FROM THE ACTUAL MAX
#include <iostream>
#include <math.h>
#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>

using namespace std;

int numTotalGraphs = 1000; //total number of graphs
int n;
int nBeginCh4 = 4300; //Beginning of cathode minimum
int nBeginCh3 = 7600; //Beginning of anode maximum
int sampleSize = 3000; //number of points to find median of graph offsets from 0
double cut = 390.;
const int batchArraySize = 6;
int batchArray [batchArraySize] = {10, 20, 25, 50, 100, 125};
double tauelec_preampA = 91.5031*1e-6; // PRE_AMP FOR ANODE
double tauelec_preampB = 43.4835*1e-6; //PRE_AMP FOR CATHODE
double gain_AoverB = 0.89512;

double getCorrection(double fittedDriftTime, double tauelec, double taulife){
  double laura = (1.-exp(-(fittedDriftTime)*(1./tauelec + 1./taulife)))/(fittedDriftTime*(1./tauelec + 1./taulife));
  return laura;
}

int maxNValue(double *yValues, int nBegin, int nPoints){
  double max = 0;
  int maxn;
  for (int i=0; i<nPoints; i++){
    if (std::fabs(yValues[i+nBegin]) > std::abs(max)){
      max = yValues[i+nBegin];
      maxn = i+nBegin;
    }
  }
  return maxn;
}

double variance(double *array, double mean, int size){
  double var = 0;
  for (int i=0; i<size; i++){
    var += (array[i]-mean)*(array[i]-mean)/(size-1);
  } 
  return var;
}

double t2Method3(double *yValues, double *xValues, int nMax)
{
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

double findMedian(double *array, int sizeOfArray) { 
    const int size = sizeOfArray;
    int index[size];
    TMath::Sort(size, array, index, false); 
    if (size % 2 != 0){
       return array[index[(size-1)/2]]; 
    }  
    else{
      return (array[index[(size-2)/2]] + array[index[size/2]])/2.0; 
    }
}

double logistic(int x, double y)
{
  if (x <= cut){
    return y;
  }
  else{
    return 0.0;
  }  
}

void fourierConv(double *original, double *fourier, int size){
  int numPoints = size;
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &numPoints, "r2c");
  fftr2c->SetPoints(original);
  fftr2c->Transform();  
  double *re = new double[numPoints];
  double *im = new double[numPoints];
  fftr2c->GetPointsComplex(re,im);
  for (int i=0; i<numPoints; i++){
    re[i] = logistic(i,re[i]);
    im[i] = logistic(i,im[i]);
  }
  TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &numPoints, "c2r");  
  fftc2r->SetPointsComplex(re,im);
  fftc2r->Transform();
  fftc2r->GetPoints(fourier);
  for (int i=0; i<numPoints; i++){
    fourier[i] = fourier[i]/numPoints;
  }
}

void leastSqrs(int nPoints, double *xValues, double *yValues, double *a, double *b, double *leastSquareYArray)
{
  double xi;
  double yi;
  double sumx = 0.0;
  double sumy = 0.0;
  double sumxy = 0.0;
  double sumx2 = 0.0;
  for (int i=0; i<nPoints; i++){
    xi = xValues[i];
    yi = yValues[i];
    if (std::isnan(yi)){
      yi = 0.0;
    }
    sumx += xi;
    sumy += yi;
    sumxy += xi*yi;
    sumx2 += xi*xi;
  }

  double denom = (nPoints*sumx2 - sumx*sumx);
  *a = (sumy*sumx2 - sumx*sumxy)/denom;
  *b = (nPoints*sumxy - sumx*sumy)/denom;
  for (int i=0; i<nPoints; i++){
    leastSquareYArray[i] = (*a) + (*b)*xValues[i];
    //cout<<"x = "<<xValues[i]<<",   y = "<<yValues[i]<<endl;
  }
  //cout<<"intersection with yaxis = "<<*a<<endl;
  //cout<<"gradient = "<<*b<<endl;
}

double t_test(double intersect, double grad, double *yValues, double *xValues, int size)
{
  //ASSUMING THAT H0: GRAD = 0 & H1: GRAD > 0
  double RSS = 0.0;
  double RSE;
  double SE;
  double xMean = 0.0;
  double xVar = 0.0;
  double tVal;
  double count = 0;
  for (int i=0; i<size; i++){
    if (!std::isnan(yValues[i])){
      RSS += pow((intersect + (grad*xValues[i]) - yValues[i]),2);
      xMean += xValues[i];
      count ++;
    }
  }
  xMean = xMean/count;
  for (int i=0; i<size; i++){
    if (!std::isnan(yValues[i])){
      xVar += pow((xValues[i]-xMean), 2.0);
    }
  }
  RSE = (RSS/(size-2));
  SE = std::sqrt(RSE/xVar);
  tVal = grad/SE;
  //cout<<"The t value is = "<<tVal<<endl;
  return tVal;
}

void totalGraphs(TFile *file, vector<double> &vec, int begin, int numIterations)
{
  double *total = new double[n]();
  for (int igraph=begin; igraph<begin+numIterations; igraph++){
    string name = Form("graph%i", igraph);
    TGraph *g = (TGraph*)file->Get(name.c_str());
    double *y = g->GetY();
    for (int ipoint = 0; ipoint<n; ipoint++){
      total[ipoint] += y[ipoint]/numIterations;
    }
  }

  double median = findMedian(total, sampleSize);
  //cout<<"median = "<<median<<endl;
  for (int i=0; i<n; i++){
    total[i] = total[i] - median;
  }

  double *fourier = new double[n]();
  fourierConv(total, fourier, n);
  vec.clear();
  for (int i=0; i<n; i++){
    vec.push_back(fourier[i]);
  }
}

double lifetime(double Qa, double Qc, double t1, double t2, double t3)
{
  double tau = (t2+(t1+t3)/2)/log(fabs(Qa/Qc));
  return tau;
}

void lifetimeArray(vector<vector<double>> &ch3VecVec, vector<vector<double>> &ch4VecVec, double *lifetimes, double *x, double *xLifetimeArray, int numIterations, double *t1Values, double *t2Values, double *t3Values, double *QaValues, double *QcValues)
{
  int maxnCh4;
  int maxnCh3;
  double Qc;
  double Qa;
  double t1;
  double t2;
  double t3;
  double tau;
  double *ch3;
  double *ch4;
  double correctionA;
  double correctionC;
  double newQa;
  double newQc;
  int count;
  
  for (int i=0; i<numIterations; i++){
    count = 0;
    ch3 = (double*)ch3VecVec[i].data();
    ch4 = (double*)ch4VecVec[i].data();
    
    maxnCh4 = maxNValue(ch4, nBeginCh4, (n-nBeginCh4));
    maxnCh3 = maxNValue(ch3, nBeginCh3, (n-nBeginCh3));
    Qc = ch4[maxnCh4];
    Qa = ch3[maxnCh3];
    Qc *= gain_AoverB;
    t1 = x[maxnCh4];
    t2 = t2Method3(ch3, x, maxnCh3)-t1;
    t3 = x[maxnCh3];
    t1Values[i] = t1;
    t2Values[i] = t2;
    t3Values[i] = t3;
    QaValues[i] = Qa;
    QcValues[i] = Qc;
    tau = lifetime(Qa, Qc, t1, t2, t3);
    while (count<20){
      correctionA = getCorrection(t3, tauelec_preampA, tau);
      correctionC = getCorrection(t1, tauelec_preampB, tau);
      newQa = Qa/correctionA;
      newQc = Qc/correctionC;
      tau = lifetime(newQa, newQc, t1, t2, t3);      
      count++; 
    }    
    tau *= 1000;
    lifetimes[i] = tau;
    xLifetimeArray[i] = (i+0.5)*(numTotalGraphs/numIterations);
    //cout<<"Qa = "<<Qa<<",    Qc = "<<Qc<<",    t1 = "<<t1<<",    t2 = "<<t2<<",    t3 = "<<t3<<",   lifetime = "<<lifetimes[i]<<endl;
  }
}


vector<double> multiLifetime(int numOfGraphsInBatch, TGraph *graphTime)
{
  int numIter = numTotalGraphs/numOfGraphsInBatch;
  TFile *ch3 = new TFile("/home/gmattimoe/data/K-890GK-800GA800A1000_50.100.200V.cm_FibreIn_23.08.ch3.traces.root", "read");
  TFile *ch4 = new TFile("/home/gmattimoe/data/K-890GK-800GA800A1000_50.100.200V.cm_FibreIn_23.08.ch4.traces.root", "read");
  TGraph *g1 = (TGraph*)ch3->Get("graph1"); //importing graph
  n = g1->GetN(); //number of points in graph
  double *x = g1->GetX(); // time values
  
  vector<vector<double>> vecVecCh3(numIter);
  vector<vector<double>> vecVecCh4(numIter);
  vector<double> yVector;
  yVector.clear();
  for (int i=0; i<numIter; i++){
    totalGraphs(ch3, yVector, (numOfGraphsInBatch*i)+1, numOfGraphsInBatch);
    vecVecCh3[i]= yVector;
  }
  yVector.clear();
  for (int i=0; i<numIter; i++){
    totalGraphs(ch4, yVector, (numOfGraphsInBatch*i)+1, numOfGraphsInBatch);
    vecVecCh4[i] = yVector;
  }  
  
  double *t1Array = new double[numIter]();
  double *t2Array = new double[numIter]();
  double *t3Array = new double[numIter]();
  double *QaArray = new double[numIter]();
  double *QcArray = new double[numIter]();
  double *logQaQc = new double[numIter]();
  double *lifetimeArr = new double[numIter]();
  double *xLifetimeArr = new double[numIter]();
  lifetimeArray(vecVecCh3, vecVecCh4, lifetimeArr, x, xLifetimeArr, numIter, t1Array, t2Array, t3Array, QaArray, QcArray);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GRAPH TIME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  graphTime->Set(numIter);
  for(int i=0; i<numIter; i++){
    graphTime->SetPoint(i, xLifetimeArr[i], lifetimeArr[i]);
  }
  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GRAPH TIME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  for (int i=0; i<numIter; i++){
    logQaQc[i] = 1/log(fabs(QaArray[i]/QcArray[i]));
  }
  
  double sum = 0.0;
  double totalNumber = 0.0;
  for (int i=0; i<numIter; i++){
    if (!std::isnan(lifetimeArr[i])){
      sum += lifetimeArr[i];
      totalNumber = totalNumber + 1.0;
    }
  }  
  double meanLifetime = sum/totalNumber;
  double a; //yIntersection
  double b; //Gradient
  double *leastSqrY = new double[numIter]();
  leastSqrs(numIter, xLifetimeArr, lifetimeArr, &a, &b, leastSqrY); //calculating the least squares values at the points in xLifetimeAarray
  
  double stdev = sqrt(variance(lifetimeArr, meanLifetime, numIter)); //calculating the standard deviation
  
  
  vector<double> values(4);
  values[0] = a;
  values[1] = b;
  values[2] = meanLifetime;
  values[3] = stdev;
  
  return values;
}

void multiBatchMultiLifetime(){
  TCanvas *c2 = new TCanvas("c1", "multipads", 1200,700);
  c2->Divide(3,2);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle("Multiple Lifetime Values for Different Batch Sizes");
  TGraph *batchOf5Graph = new TGraph();
  TGraph *batchOf10Graph = new TGraph();
  TGraph *batchOf20Graph = new TGraph();
  TGraph *batchOf25Graph = new TGraph();
  TGraph *batchOf50Graph = new TGraph();
  TGraph *batchOf100Graph = new TGraph();
  TGraph *batchOf125Graph = new TGraph();
  TGraph *batchOf200Graph = new TGraph();
  vector<TGraph*> graphVector{batchOf10Graph, batchOf20Graph, batchOf25Graph, batchOf50Graph, batchOf100Graph, batchOf125Graph};
  double intersects [batchArraySize];
  double gradients [batchArraySize];
  double averageLifetimes [batchArraySize];
  double stdev [batchArraySize];
  
  vector<double> valueVector(4);
  for (int i=0; i<batchArraySize; i++){
    valueVector.clear();
    valueVector = multiLifetime(batchArray[i], graphVector[i]);
    intersects[i] = valueVector[0];
    gradients[i] = valueVector[1];
    averageLifetimes[i] = valueVector[2];
    stdev[i] = valueVector[3]; 
  }
  
  for (int i=0; i<batchArraySize; i++){
    cout<<"The batch size = "<<batchArray[i]<<",   intersect = "<<intersects[i]<<",   gradient = "<<gradients[i]<<",   average lifetime = "<<averageLifetimes[i]<<",   standard deviation = "<<stdev[i]<<endl;
  }
  

  int colorArray[8] = {46, 38, 8, 1, 4, 2, 12, 32};
  for (int i=0; i<batchArraySize; i++){
    graphVector[i]->SetMarkerStyle(8);
    graphVector[i]->SetMarkerSize(0.8);
    graphVector[i]->SetLineStyle(2);
    graphVector[i]->SetLineWidth(1);
    graphVector[i]->SetMarkerColor(colorArray[i]);
    graphVector[i]->SetLineColor(colorArray[i]);
    c2->cd(i);
    gPad->Update();
    graphVector[i]->Draw("ALP");
    if (i==0 || i==1){
      mg1->Add(graphVector[i], "P");
      graphVector[i]->SetMarkerStyle(23);
    }
    else{
      mg1->Add(graphVector[i], "PL");
      graphVector[i]->SetMarkerStyle(8);
    }
  }
  
  cout<<"and wE got This FAR "<<endl;
  mg1->Draw("ALP");
  cout<<"and we have even gone This FAR"<<endl;
  mg1->GetXaxis()->SetTitle("Middle Waveform Number");
  mg1->GetYaxis()->SetTitle("Lifetime ms");
  //mg1->GetXaxis()->SetLimits(-0.00005,0.00055);
  //mg1->GetYaxis()->SetRangeUser(0.0, 0.001);
  gPad->SetGrid(1, 1); 
  gPad->Update(); 
  
  auto legend = new TLegend(0.1,0.75,0.45,0.9);
  legend->SetNColumns(3);
  legend->SetTextSize(0.025);
  //legend->AddEntry(batchOf5Graph, "Batch of 5", "LP");
  legend->AddEntry(batchOf10Graph, "Batch of 10", "P");
  legend->AddEntry(batchOf20Graph, "Batch of 20", "P");
  legend->AddEntry(batchOf25Graph, "Batch of 25", "LP");
  legend->AddEntry(batchOf50Graph, "Batch of 50", "LP");
  legend->AddEntry(batchOf100Graph, "Batch of 100", "LP");
  legend->AddEntry(batchOf125Graph, "Batch of 125", "LP");
  //legend->AddEntry(batchOf200Graph, "Batch of 200", "LP");
  legend->Draw(); 
  
   
  
}