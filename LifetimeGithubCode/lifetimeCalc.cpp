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
double cut = 390.;
double tauelec_preampA = 91.5031*1e-6; // PRE_AMP FOR CATHODE (ASSUMING AMPS ARE "SWITCHED")
double tauelec_preampB = 43.4835*1e-6; //PRE_AMP FOR ANODE (ASSUMING AMPS ARE "SWITCHED")
double gain_AoverB = 0.89512;
int sampleSize = 3000;  
int nBeginCh4 = 4400; //Beginning of cathode minimum
int nBeginCh3 = 7800; //Beginning of anode maximum


double getCorrection(double fittedDriftTime, double tauelec, double taulife){
  double laura = (1.-exp(-(fittedDriftTime)*(1./tauelec + 1./taulife)))/(fittedDriftTime*(1./tauelec + 1./taulife));
  return laura;
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


void lifetimeCalc(){
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
  double *fourierCh3 = new double[n]();
  double *fourierCh4 = new double[n]();
  
  fourierConv(totalCh3, fourierCh3, n);
  fourierConv(totalCh4, fourierCh4, n);
  int maxnCh4 = maxNValue(fourierCh4, nBeginCh4, (n-nBeginCh4));
  int maxnCh3 = maxNValue(fourierCh3, nBeginCh3, (n-nBeginCh3));
  double Qc = fourierCh4[maxnCh4];
  double Qa = fourierCh3[maxnCh3];
  cout<<"Qc = "<<Qc<<endl;
  Qc *= gain_AoverB;
  int count = 0;
  double newQa;
  double newQc;
  double correctionA;
  double correctionC;
  double t1 = x[maxnCh4];
  double t2 = t2Method3(fourierCh3, x, maxnCh3)-t1;
  double t3 = x[maxnCh3];
  double tau = lifetime(Qa, Qc, t1, t2, t3); //initial estimate of lifetime
  while (count<20){
    correctionA = getCorrection(t3, tauelec_preampA, tau);
    correctionC = getCorrection(t1, tauelec_preampB, tau);
    newQa = Qa/correctionA;
    newQc = Qc/correctionC;
    tau = lifetime(newQa, newQc, t1, t2, t3);      
    count++; 
  }
  cout<<"lifetime = "<<tau<<endl;
  
  
  TCanvas *c2 = new TCanvas("c1", "multipads", 900,700);
  TGraph *anodeGraph = new TGraph(n, x, fourierCh3);
  TGraph *cathodeGraph = new TGraph(n, x, fourierCh4);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(anodeGraph);
  mg1->Add(cathodeGraph);
  mg1->Draw("ALP");
}