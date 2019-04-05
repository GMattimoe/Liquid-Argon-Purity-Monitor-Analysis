//adaptation of cathode plot to get times withough using filename times

#include<sys/types.h>
#include<dirent.h>
#include<unistd.h>
#include "functions.h"

int lengthOf100 = 20;
double numSecPerDay = 86400;

int findAllTimes(string folder, string times[3000], int minutes[3000]);

double timeToDouble(string timeAndDate);

void removingGapsInTime(map<string, double> map,  std::map<string, double> &mapPointer);

map<string, double> mapAllTimes(string folder);

string readTime(string folder, string filename);

//void plotCathodeVsTime(string basename="2018September11Vacumm/NewLampPosition/Gold/ContinuousExposition/"){
void newCathodeTimePlot(string basename = "/unix/dune/purity/2019February21Vacuum/") {
  string title = "Cathode data";
  string material[] = {"Gold","Silver"};
  string titleout = ("Vacuum_Cathode_Amplitude_With_Gap_" + material[0] +"_"+ material[1]).c_str();  
  string times[3000];
  int minutes[3000];
  double x[2][3000];
  double xmin;  
  double yPeak[2][3000];
  double timePeak[2][3000];
  //double yFitPeak[2][3000];
  TGraph * gPeak[2];
  TGraph * tPeak[2];
  static const char* funcNames[] = {"GoldFunc","SilverFunc"};//{"Silver_exp", "Silver_1/x", "Gold_exp", "Gold_1/x"};
  TF1 *funcs[2];
  //TGraph * gFitPeak[2];
  TCanvas * c1 = new TCanvas("c");
  TCanvas * c2 = new TCanvas("c2");
  TLegend * leg = new TLegend(0.6, 0.11, 0.9, 0.25);
  double tauelecK = 43.4835 * 1e-6;
  vector<map<string,double>> vecMap(2);
  
  for (int imat = 0; imat < 2; imat++) {
    int countpoints = 0;
    string folder = basename + material[imat];
    cout << folder << endl;
    cout<<"At line "<<__LINE__<<endl;
    vecMap[imat] = mapAllTimes(folder);
    removingGapsInTime(vecMap[imat], vecMap[imat]);
    cout<<"At line "<<__LINE__<<endl;
    //map<string, string> silverMap;
    //mapAllTimes(folder, silverMap);
    std::map<string,double>::iterator iter;
    for (iter = vecMap[imat].begin(); iter!=vecMap[imat].end(); ++iter) {
      cout<<"At line "<<__LINE__<<endl;
      string filename = iter->first;
      cout << filename << endl;
      TFile *fin = new TFile((folder + "/" + filename).c_str(), "read");
      double waveTimes = iter->second;
      TGraph * gtemp = (TGraph * ) fin->Get("zeroedAvgSmooth;1");

      //TF1 * func = new TF1("func", fittingFunction, 0., 0.9E-3, 4);
      //func->SetParameters(5, 90, -10, 1);
      //func->SetParLimits(1, tauelecK * 1E6 * 0.95, tauelecK * 1E6 * 1.05);
      //func->SetParName(0, "#sigma");
      //func->SetParName(1, "#tau_{D}");
      //func->SetParName(2, "a");
      //func->SetParName(3, "t_0");
      //gtemp->Fit("func", "QC", "RQN0", -0.1E-3, 0.5E-3);
      //yFitPeak[imat][countpoints] = func->GetMinimum() * (-1);
      yPeak[imat][countpoints] = TMath::MinElement(gtemp->GetN(), gtemp->GetY()) * (-1);
      timePeak[imat][countpoints] = gtemp->GetX()[TMath::LocMin(gtemp->GetN(), gtemp->GetY())] * (1e6);
      x[imat][countpoints] = waveTimes;
      countpoints++;
      //delete func;
      delete gtemp;
      fin->Close();
      delete fin;
    }
    
    gPeak[imat] = new TGraph(countpoints, x[imat], yPeak[imat]);
    xmin = TMath::MinElement(gPeak[imat]->GetN(), gPeak[imat]->GetX());
    for(int i=0; i<countpoints; i++){
      x[imat][i] -= xmin;
    }
    gPeak[imat] = new TGraph(countpoints, x[imat], yPeak[imat]);
    gPeak[imat]->SetLineColor(imat + 1);
    gPeak[imat]->SetMarkerColor(imat + 1);
    gPeak[imat]->SetMarkerStyle(2);
    gPeak[imat]->SetLineWidth(2);
    //leg->AddEntry(gPeak[imat], material[imat].c_str(), "l");
    gPeak[imat]->SetTitle((title + ";Time [s];Peak [V]").c_str());
    /*
    funcs[imat] = new TF1(funcNames[imat],"[0]+[1]/(x+[2])");
    gPeak[imat]->Fit(funcs[imat], "M");
    gPeak[imat]->GetFunction(funcNames[imat])->SetLineColor(imat+1);
    */
    /*
    funcs[imat*2] = new TF1(funcNames[imat*2],"[0]+exp([1]+[2]*x)");
    gPeak[imat*2]->Fit(funcs[imat*2], "M");
    gPeak[imat*2]->GetFunction(funcNames[imat*2])->SetLineColor(imat+1);
    funcs[imat*2 + 1] = new TF1(funcNames[imat*2+1],"[0]+([1]/x)");
    gPeak[imat*2 + 1]->Fit(funcs[imat*2+1], "M");
    gPeak[imat*2 + 1]->GetFunction(funcNames[imat*2+1])->SetLineColor(imat+1);
    */    
    //    g[imat]->GetYaxis()->SetRangeUser(0, 4);
/*
    gFitPeak[imat] = new TGraph(countpoints, x[imat], yFitPeak[imat]);
    gFitPeak[imat]->SetLineColor(imat + 1);
    gFitPeak[imat]->SetLineWidth(2);
    gFitPeak[imat]->SetMarkerColor(imat + 1);
    gFitPeak[imat]->SetMarkerStyle(2);
    gFitPeak[imat]->SetTitle((title + ";Time [s];Fitted Peak [V]").c_str());
*/
    tPeak[imat] = new TGraph(countpoints, x[imat], timePeak[imat]);
    tPeak[imat]->SetLineColor(imat + 1);
    tPeak[imat]->SetMarkerColor(imat + 1);
    tPeak[imat]->SetMarkerStyle(2);
    tPeak[imat]->SetLineWidth(2);
    //leg->AddEntry(gPeak[imat], material[imat].c_str(), "l");
    tPeak[imat]->SetTitle((title + ";Time [s];Time of Peak [#mus]").c_str());    
  }



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GETS TO HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << "Out here" << endl;

  int nG, nS;
  TFile * fout = new TFile((titleout + ".root").c_str(), "recreate");
  if (gPeak[0] && gPeak[1]) {
    cout << "DONE?" << endl;
    gStyle->SetOptStat(0);
    c1->cd();
    TMultiGraph *mg1 = new TMultiGraph();
    mg1->SetTitle("Cathode Amplitude in Vacuum For Silver and Gold;Time [Days];Peak [V]");
    mg1->Add(gPeak[0], "P");
    mg1->Add(gPeak[1], "P");
    mg1->Draw("ALP");
    auto legend1 = new TLegend(0.65,0.75,0.9,0.9);
    //legend->SetNColumns(3);
    //legend->SetTextSize(0.025);
    legend1->AddEntry(gPeak[1], "Silver", "P");
    legend1->AddEntry(gPeak[0], "Gold", "P");
    legend1->Draw();   
    c1->Print((titleout + ".png").c_str());
    c1->Print((titleout + ".pdf").c_str());
    c2->cd();
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->SetTitle("Cathode Time of Peak Amplitude in Vacuum For Silver and Gold;Time [Days];Time of Peak [#mus]");
    mg2->Add(tPeak[0], "P");
    mg2->Add(tPeak[1], "P");
    mg2->Draw("ALP");
    auto legend2 = new TLegend(0.65,0.45,0.9,0.6);
    //legend->SetNColumns(3);
    //legend->SetTextSize(0.025);
    legend2->AddEntry(tPeak[1], "Silver", "P");
    legend2->AddEntry(tPeak[0], "Gold", "P");
    legend2->Draw();   
    c2->Print((titleout + "peaktime.png").c_str());
    c2->Print((titleout + "peaktime.pdf").c_str());
    

  }
  
  cout<<"At line "<<__LINE__<<endl;
  if (gPeak[0]) {
    //gPeak[0]->Draw("AP");
    //gFitPeak[0]->Draw("AP");
    gPeak[0]->Write("gPeak_gold");
    tPeak[0]->Write("timePeak_gold");
    //gFitPeak[0]->Write("gFitPeak_gold");
  }
  
  if (gPeak[1]) {
    gPeak[1]->Write("gPeak_silver");
    tPeak[1]->Write("timePeak_silver");
    //gFitPeak[1]->Write("gFitPeak_silver");
  }
  

}

int findNumTimes(string folder) {
  DIR * dp;
  dirent * d;
  int count = 0;
  string pattern1 = ".ch4.traces_averages.root";
  if ((dp = opendir(folder.c_str())) == NULL)
    perror("opendir");

  while ((d = readdir(dp))) {
      count++;
  }
  return count;

}

string readTime(string folder, string filename){ 
  //reads the time and date from the .times file returns a string in the form - "year_month_day_hour_minute_second"
  string locAndName = folder + "/" + filename;
  ifstream inFile(locAndName);
  std::string line;
  std::getline(inFile, line);
  string tokens[10];
  std::string::size_type pos;
  int count=0;
  string delimiter = ", ";
  line.erase(0,1);
  while((pos = line.find(delimiter)) != std::string::npos){
    tokens[count] = line.substr(0,pos);
    line.erase(0,pos+delimiter.length());
    count++;
  }
  string dateAndTime = tokens[0];
  for (int i=1; i<count; i++){
    dateAndTime = tokens[i] + "_" + dateAndTime;
  }
  return dateAndTime;
}

double timeToDouble(string timeAndDate){
  std::string::size_type pos;
  struct tm tm;
  string tokens [6];
  int itoken = 0;
  while((pos = timeAndDate.find("_")) != std::string::npos){
    tokens[itoken] = timeAndDate.substr(0,pos);
    timeAndDate.erase(0,pos+1);
    itoken++;
  }
  tokens[itoken] = timeAndDate;
  tm.tm_sec = atof(tokens[5].c_str());
  tm.tm_min = atoi(tokens[4].c_str());
  tm.tm_hour = atoi(tokens[3].c_str());
  tm.tm_mday = atoi(tokens[2].c_str());
  tm.tm_mon = atoi(tokens[1].c_str())-1;
  tm.tm_year = atoi(tokens[0].c_str())-1949;   
  double waveTimes = mktime( & tm);  
  return waveTimes/numSecPerDay;
}

map<string, double> mapAllTimes(string folder) {
  map<string, double> map;
  double doubleTime;
  //MAP THE  filenames to the time string of the time
  DIR * dp;
  dirent * d;
  string pattern1 = ".ch4.traces_averages.root";
  string pattern2 = ".ch4.traces.times";
  string timeAndDate;
  string timeFileName;
  //cout<<folder<<endl;  
  if ((dp = opendir(folder.c_str())) == NULL){
    perror("opendir");
  }
  while ((d = readdir(dp)) != NULL) {
    if (!strcmp(d->d_name, ".") || !strcmp(d->d_name, ".."))
      continue;
    std::string temps = d->d_name;
    if (temps.find(pattern1) != std::string::npos) {
      std::string::size_type i = temps.find(pattern1);
      timeFileName = temps.substr(0, i);
      timeFileName += pattern2;
      timeAndDate = readTime(folder, timeFileName);
      doubleTime = timeToDouble(timeAndDate);
      map.insert(std::pair<string,double>(temps, doubleTime));
    }
  }
  return map;
}

void removingGapsInTime(map<string, double> map, std::map<string, double> &mapPointer){
  std::map<string,double>::iterator iter;
  string pattern = "Second_Set";
  double maxTime = 0.0;
  double minTime = 1e10;
  for (iter = map.begin(); iter!=map.end(); ++iter) {
    if(iter->first.find(pattern) != std::string::npos){
      if(iter->second<minTime)
      minTime = iter->second;
    }
    else
      if(iter->second>maxTime){
        maxTime = iter->second;
      }
  }
  cout<<"the minimum time = "<<minTime<<endl;
  double newTime;
  for (iter = map.begin(); iter!=map.end(); ++iter) {
    if(iter->first.find(pattern) != std::string::npos){
      newTime = iter->second - minTime + maxTime;
      //mapPointer.find(iter->first)->second = newTime;
    }
  }
}

