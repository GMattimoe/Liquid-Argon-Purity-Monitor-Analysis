//SPECTRUM ANALYSIS TO FIND PEAKS IN OSCILLATIONS
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>
using namespace std;

int numGraphs = 1000; //number of graphs to be averaged
int n;
int sampleSize = 3000;
int subFourierSize = 4000;
double numStDev = 1.0;
int beginFit1 = 100;
int endFit1 = 300;
int beginFit2 = 300;
int endFit2 = 500;
int beginFit3 = 500;
int endFit3 = 720;
int beginFit4 = 720;
int endFit4 = 900;
int bufferLength = 10;


void fourierTrans(double *original, double *spectrum, int inputSize, int outputSize){
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &inputSize, "r2c");
  fftr2c->SetPoints(original);
  fftr2c->Transform();  
  double *re = new double[outputSize];
  double *im = new double[outputSize];
  fftr2c->GetPointsComplex(re,im);
  for (int i=0; i<outputSize; i++){
    spectrum[i] = std::sqrt(((re[i]*re[i]) + (im[i]*im[i])));
  }
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

double standardDeviation(double *array, double mean, int size){
  double var = 0;
  for (int i=0; i<size; i++){
    var += (array[i]-mean)*(array[i]-mean)/(size-1);
  }
  double sd = sqrt(var);  
  return sd;
}

void spectralSmoothing(double *signalSpectrum, double *noiseSpectrum, double *smoothedSignalSpectrum, int signalArraySize, int noiseArraySize, double noiseTolerance){
  int noiseInt;
  int noiseCount;
  int signalCount;
  double newSpectrumValue;
  double startValue;
  double endValue;
  for (int i=0; i<signalArraySize; i++){
    noiseCount = 0;
    signalCount = 0;
    noiseInt = i*noiseArraySize/signalArraySize;
    if (i<100){
      smoothedSignalSpectrum[i] = signalSpectrum[i];
    }
    else{
      while (noiseSpectrum[noiseInt] > noiseTolerance){
        noiseCount++;
        noiseInt++;
      }
      if (noiseCount>0){
        signalCount = noiseCount*signalArraySize/noiseArraySize;
        for (int j=0; j<signalCount; j++){
          newSpectrumValue = signalSpectrum[i] + (signalSpectrum[i+signalCount]-signalSpectrum[i])*j/signalCount;
          smoothedSignalSpectrum[i+j] = newSpectrumValue;
        } 
        i += (signalCount-1);
      }
      else{
        smoothedSignalSpectrum[i] = signalSpectrum[i];
      }
    }  
  }
}

void logSpectralRemoval(double *signalSpectrum, double *noiseSpectrum, TGraph *removedSpectralSignal, int signalArraySize, int noiseArraySize, double noiseTolerance, double timePeriod){
  int noiseInt;
  int noiseCount;
  int pointInt;
  int signalCount;
  double newSpectrumValue;
  double startValue;
  double endValue;
  double frequency;  
  for (int i=0; i<signalArraySize; i++){
    pointInt = removedSpectralSignal->GetN();
    noiseCount = 0;
    signalCount = 0;
    noiseInt = i*noiseArraySize/signalArraySize;
    //if (i<100){
    //  removedSpectralSignal->SetPoint(pointInt, i, log(signalSpectrum[i]+1));
    //}
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
        frequency = (double)i/timePeriod;      
        removedSpectralSignal->SetPoint(pointInt, frequency, log(signalSpectrum[i]+1));
      }
    }  
  }
}


void spectrumAnalysis(){
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
  double *subTotalCh3 = new double[subFourierSize]();
  double *subTotalCh4 = new double[subFourierSize]();
  for (int i=0; i<subFourierSize; i++){
    subTotalCh3[i] = totalCh3[i+(n-1)-subFourierSize];
    subTotalCh4[i] = totalCh4[i+(n-1)-subFourierSize];
  }  
  int specSize = n/2 + 1;
  int specSubSize = subFourierSize/2 + 1;
  double *spectrumCh3 = new double[specSize]();
  double *spectrumCh4 = new double[specSize]();
  double *spectrumSubCh3 = new double[specSubSize]();
  double *spectrumSubCh4 = new double[specSubSize]();
  fourierTrans(totalCh3, spectrumCh3, n, specSize);
  fourierTrans(totalCh4, spectrumCh4, n, specSize);
  fourierTrans(subTotalCh3, spectrumSubCh3, subFourierSize, specSubSize);
  fourierTrans(subTotalCh4, spectrumSubCh4, subFourierSize, specSubSize);
  double timePeriod = (x[n-1] - x[0])*1000;
  double timePeriodSub = (x[subFourierSize-1] - x[0])*1000;
  double *frequency = new double[specSize]();
  double *frequencySub = new double[specSubSize]();
  for (int i=0; i<specSize; i++){
    frequency[i] = (double)i/timePeriod;
  }
  for (int i=0; i<specSubSize; i++){
    frequencySub[i] = (double)i/timePeriodSub;
  }
  
  
  
  double noiseMeanCh3 = 0.;
  for (int i=0; i<specSubSize; i++){
    noiseMeanCh3 += spectrumSubCh3[i]/specSubSize;
  }
  double noiseMeanCh4 = 0.;
  for (int i=0; i<specSubSize; i++){
    noiseMeanCh4 += spectrumSubCh4[i]/specSubSize;
  }  
  
  double noiseStandardDeviationCh3 = standardDeviation(spectrumSubCh3, noiseMeanCh3, specSubSize);
  double noiseStandardDeviationCh4 = standardDeviation(spectrumSubCh4, noiseMeanCh4, specSubSize);
  cout<<noiseStandardDeviationCh3<<endl;
  double noiseToleranceCh3 = noiseMeanCh3 + numStDev*noiseStandardDeviationCh3;
  double noiseToleranceCh4 = noiseMeanCh4 + numStDev*noiseStandardDeviationCh4;
  double *smoothedSpectrumCh3 = new double[specSize]();
  spectralSmoothing(spectrumCh3, spectrumSubCh3, smoothedSpectrumCh3, specSize, specSubSize, noiseStandardDeviationCh3);
  
  TGraph *logSpectralRemovalGraphCh3 = new TGraph();
  cout<<logSpectralRemovalGraphCh3->GetN()<<endl;
  logSpectralRemoval(spectrumCh3, spectrumSubCh3, logSpectralRemovalGraphCh3, specSize, specSubSize, noiseToleranceCh3, timePeriod);

  TGraph *logSpectralRemovalGraphCh4 = new TGraph();
  cout<<logSpectralRemovalGraphCh4->GetN()<<endl;
  logSpectralRemoval(spectrumCh4, spectrumSubCh4, logSpectralRemovalGraphCh4, specSize, specSubSize, noiseToleranceCh4, timePeriod);
  
  //cout<<specSize<<endl;
  //cout<<logSpectralRemovalGraphCh3->GetN()<<endl;
  double *logSpectrumCh3 = new double[specSize]();
  double *logSpectrumCh4 = new double[specSize]();
  double *logSmoothedSpectrumCh3 = new double[specSize]();
  for (int i=0; i<specSize; i++){
    logSpectrumCh3[i] = log(spectrumCh3[i]+1);
    logSpectrumCh4[i] = log(spectrumCh4[i]+1);
    logSmoothedSpectrumCh3[i] = log(smoothedSpectrumCh3[i]+1);
  }
  double *logSpectrumSubCh3 = new double[specSubSize]();
  double *logSpectrumSubCh4 = new double[specSubSize]();
  for (int i=0; i<specSubSize; i++){
    logSpectrumSubCh3[i] = log(spectrumSubCh3[i]+1);
    logSpectrumSubCh4[i] = log(spectrumSubCh4[i]+1);
  }
  

  
  
  
  
  TCanvas *c2 = new TCanvas("c1", "multipads", 900,700);
  c2->Divide(1,2);
  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();
  TGraph *fullSpectrumGraphCh3 = new TGraph(specSize, frequency, spectrumCh3);
  TGraph *logSpectrumGraphCh3 = new TGraph(specSize, frequency, logSpectrumCh3);
  TGraph *fullSpectrumGraphCh4 = new TGraph(specSize, frequency, spectrumCh4);
  TGraph *logSpectrumGraphCh4 = new TGraph(specSize, frequency, logSpectrumCh4); 
  TGraph *logSmoothedSpectrumGraphCh3 = new TGraph(specSize, frequency, logSmoothedSpectrumCh3);
  
  fullSpectrumGraphCh3->SetMarkerColor(9);
  fullSpectrumGraphCh4->SetMarkerColor(9);
  logSpectrumGraphCh3->SetMarkerColor(9);
  logSpectrumGraphCh4->SetMarkerColor(9);
  
  fullSpectrumGraphCh3->SetLineColor(9);
  fullSpectrumGraphCh4->SetLineColor(9);
  logSpectrumGraphCh3->SetLineColor(9);
  logSpectrumGraphCh4->SetLineColor(9);
  
  fullSpectrumGraphCh3->SetTitle("Full Frequency Spectrum of Anode Signal");
  fullSpectrumGraphCh4->SetTitle("Full Frequency Spectrum of Cathode Signal");
  logSpectrumGraphCh3->SetTitle("Log Frequency Spectrum of Anode Signal");
  logSpectrumGraphCh4->SetTitle("Log Frequency Spectrum of Catode Signal");
  

  TGraph *logSpectrumSubGraphCh3 = new TGraph(specSubSize, frequencySub, logSpectrumSubCh3);
  TGraph *logSpectrumSubGraphCh4 = new TGraph(specSubSize, frequencySub, logSpectrumSubCh4);  
  
  logSpectrumSubGraphCh3->SetMarkerColor(1);
  logSpectrumSubGraphCh4->SetMarkerColor(1);
  
  logSpectrumSubGraphCh3->SetLineColor(1);
  logSpectrumSubGraphCh4->SetLineColor(1);
  
  logSpectrumSubGraphCh3->SetTitle("Log Frequency Spectrum of Sample of Anode Signal");
  logSpectrumSubGraphCh4->SetTitle("Log Frequency Spectrum of Sample of Catode Signal");      
  

  /*
  c2->cd(1);
  logSpectrumGraphCh3->Draw();

  c2->cd(3);
  logSpectrumGraphCh4->Draw(); 
  
  c2->cd(2);
  logSpectrumSubGraphCh3->Draw();
  
  c2->cd(4);
  logSpectrumSubGraphCh4->Draw();   
  */
  

  TF1 *fit1Ch3 = new TF1("fit1", "pol6", beginFit1/timePeriod, (endFit1+bufferLength)/timePeriod);
  TF1 *fit2Ch3 = new TF1("fit2", "pol6", beginFit2/timePeriod, (endFit2+bufferLength)/timePeriod);
  TF1 *fit3Ch3 = new TF1("fit3", "pol6", beginFit3/timePeriod, (endFit3+bufferLength)/timePeriod);
  TF1 *fit4Ch3 = new TF1("fit4", "pol6", beginFit4/timePeriod, (endFit4+bufferLength)/timePeriod);
  logSpectralRemovalGraphCh3->Fit(fit1Ch3, "CQR");
  logSpectralRemovalGraphCh3->Fit(fit2Ch3, "CQR+");
  logSpectralRemovalGraphCh3->Fit(fit3Ch3, "CQR+");
  logSpectralRemovalGraphCh3->Fit(fit4Ch3, "CQMR+");

  TF1 *fit1Ch4 = new TF1("fit1", "pol6", beginFit1/timePeriod, (endFit1+bufferLength)/timePeriod);
  TF1 *fit2Ch4 = new TF1("fit2", "pol6", beginFit2/timePeriod, (endFit2+bufferLength)/timePeriod);
  TF1 *fit3Ch4 = new TF1("fit3", "pol6", beginFit3/timePeriod, (endFit3+bufferLength)/timePeriod);
  TF1 *fit4Ch4 = new TF1("fit4", "pol6", beginFit4/timePeriod, (endFit4+bufferLength)/timePeriod);
  logSpectralRemovalGraphCh4->Fit(fit1Ch4, "CQR");
  logSpectralRemovalGraphCh4->Fit(fit2Ch4, "CQR+");
  logSpectralRemovalGraphCh4->Fit(fit3Ch4, "CQR+");
  logSpectralRemovalGraphCh4->Fit(fit4Ch4, "CQMR+");  
  
  c2->cd(1);
  mg1->SetTitle("Log Amplitudes in the Frequency Domain of the Anode Signal With Noise Spectrum");
  mg1->Add(logSpectrumGraphCh3, "L");
  mg1->Add(logSpectrumSubGraphCh3, "L");
  //mg1->Add(logSmoothedSpectrumGraphCh3);
  mg1->Add(logSpectralRemovalGraphCh3, "P");
  //mg1->Add(fullSpectrumGraphCh3);
  mg1->Draw("ALP");
  mg1->GetXaxis()->SetLimits(0.,500);
  mg1->GetXaxis()->SetTitle("Frequency kHz");
  mg1->GetYaxis()->SetTitle("Log Amplitude");
  auto legend1 = new TLegend(0.7,0.75,0.9,0.9);
  legend1->SetTextSize(0.035);
  legend1->AddEntry(logSpectrumGraphCh3, "Anode Signal", "l");
  legend1->AddEntry(logSpectrumSubGraphCh3, "Noise Signal", "l");
  legend1->AddEntry(logSpectralRemovalGraphCh3->GetFunction("fit1"), "Fitted Line", "l");
  legend1->Draw(); 
  
  c2->cd(2);
  mg2->SetTitle("Log Amplitudes in the Frequency Domain of the Cathode Signal With Noise Spectrum");
  mg2->Add(logSpectrumGraphCh4, "L");
  mg2->Add(logSpectrumSubGraphCh4, "L");
  mg2->Add(logSpectralRemovalGraphCh4, "P");
  mg2->Draw("ALP");
  mg2->GetXaxis()->SetLimits(0.,500);
  mg2->GetXaxis()->SetTitle("Frequency kHz");
  mg2->GetYaxis()->SetTitle("Log Amplitude");
  auto legend2 = new TLegend(0.7,0.75,0.9,0.9);
  legend2->SetTextSize(0.035);
  legend2->AddEntry(logSpectrumGraphCh4, "Cathode Signal", "l");
  legend2->AddEntry(logSpectrumSubGraphCh4, "Noise Signal", "l");
  legend2->AddEntry(logSpectralRemovalGraphCh4->GetFunction("fit1"), "Fitted Line", "l");
  legend2->Draw();  
  
  
  
  
  
}