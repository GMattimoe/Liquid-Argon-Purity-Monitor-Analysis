//function to fit various model functions on cathode data

void fittingTime(){
  TFile *fin = new TFile("/home/gmattimoe/CathodeLongExposure/purityMonitorAnalysis/Vacuum_Cathode_Amplitude_With_Gap_Gold_Silver.root", "read");
  TGraph *g1 = (TGraph*)fin->Get("gPeak_gold;1");
  TGraph *s1 = (TGraph*)fin->Get("gPeak_silver;1");
  double *gX = g1->GetX();
  double *sX = s1->GetX();
  int N = g1->GetN();
  int num;
  for (int i=0; i<N; i++){
    if (gX[i]<5){
      num = i;
    }
  }
  num++;
  double *goldExcludedX = new double[num];
  double *goldExcludedY = new double[num];
  for (int i=0; i<num; i++){
    goldExcludedX[i] = gX[i];
    goldExcludedY[i] = g1->GetY()[i];
  } 
  TGraph *fullGoldGraph = new TGraph(N, gX, g1->GetY());
  TGraph *goldGraph = new TGraph(num, goldExcludedX, g1->GetY());
  TGraph *silverGraph = new TGraph(N, sX, s1->GetY());
  fullGoldGraph->SetMarkerStyle(22);
  fullGoldGraph->SetMarkerSize(0.85);
  fullGoldGraph->SetMarkerColor(1);
  
  
  silverGraph->SetMarkerStyle(7);
  silverGraph->SetMarkerSize(0.85);
  silverGraph->SetMarkerColor(1);  
  TCanvas *c2 = new TCanvas("c1", "multipads", 1200,700);
  TMultiGraph *mg = new TMultiGraph();
  
  TF1 *gold_combo = new TF1("gold_combo","[0] + [1]/([2]+x) + [3]/(x*x + [4]*x + [5])",0,5);
  goldGraph->Fit("gold_combo", "RM+","",0,400);
  goldGraph->GetFunction("gold_combo")->SetLineColor(2); 
  TF1 *full_gold_combo = new TF1("full_gold_combo","[0] + [1]/([2]+x) + [3]/(x*x + [4]*x + [5])",0,5);
  fullGoldGraph->Fit("full_gold_combo", "RM+","",0,400);
  fullGoldGraph->GetFunction("full_gold_combo")->SetLineColor(9);
  
  //TF1 *silver_combo = new TF1("silver_combo", "[0] + [1]/([2]+x) + [3]/(x*x + [4]*x + [5])");
  //TFitResultPtr bestFit = silverGraph->Fit("silver_combo", "SM+");//,"",0,400);
  //silverGraph->GetFunction("silver_combo")->SetLineColor(9);
  /*
  //EXP
  TF1 *gold_exp = new TF1("gold_exp","[0]+exp([1]+[2]*x)",0,5);
  goldGraph->Fit("gold_exp", "RM+","",0,400);
  goldGraph->GetFunction("gold_exp")->SetLineColor(46);
  TF1 *silver_exp = new TF1("silver_exp","[0]+exp([1]+[2]*x)");
  silverGraph->Fit("silver_exp", "M+");
  silverGraph->GetFunction("silver_exp")->SetLineColor(46);
    
  // 1/x
  TF1 *gold_inv = new TF1("gold_inv","[0]+[1]/([2]+x)",0,5);
  goldGraph->Fit("gold_inv", "RM+","",0,400);
  goldGraph->GetFunction("gold_inv")->SetLineColor(8);
  TF1 *silver_inv = new TF1("silver_inv","[0]+[1]/([2]+x)");
  silverGraph->Fit("silver_inv", "M+");
  silverGraph->GetFunction("silver_inv")->SetLineColor(8);
  */

/*
  //gaus
  TF1 *gold_gaus = new TF1("gold_gaus","[0]+[1]*exp(-0.5*((x-[2])/[3])^2)");
  gold_gaus->SetParameter(2,-1);
  gold_gaus->SetParameter(3,2);
  goldGraph->Fit("gold_gaus", "M+");
  goldGraph->GetFunction("gold_gaus")->SetLineColor(8); 
  TF1 *silver_gaus = new TF1("silver_gaus","[0]+[1]*exp(-0.5*((x-[2])/[3])^2)");
  silver_gaus->SetParameter(2,-1);
  silver_gaus->SetParameter(3,2);
  silverGraph->Fit("silver_gaus", "M+");
  silverGraph->GetFunction("silver_gaus")->SetLineColor(8); 
  
  TF1 *gold_combo = new TF1("gold_combo","[0]+[1]/([2]+pow(x+[3],2))");
  goldGraph->Fit("gold_combo", "M+");
  goldGraph->GetFunction("gold_combo")->SetLineColor(9); 
  
  TF1 *silver_combo = new TF1("silver_combo","[0]+[1]/([2]+pow(x+[3],2))");
  //silverGraph->Fit("silver_combo", "M+");
  //silverGraph->GetFunction("silver_combo")->SetLineColor(9);   
*/
/*
  cout<<"The best fit is with the sum of reciprocals: "<<endl;
  cout<<"p0 = "<<bestFit->Parameter(0)<<endl;
  cout<<"p1 = "<<bestFit->Parameter(1)<<endl;
  cout<<"p2 = "<<bestFit->Parameter(2)<<endl;
  cout<<"p3 = "<<bestFit->Parameter(3)<<endl;
  cout<<"p4 = "<<bestFit->Parameter(4)<<endl;
  cout<<"p5 = "<<bestFit->Parameter(5)<<endl;
  cout<<"with chi squared = "<<bestFit->Chi2()<<endl;
*/ 
  

  mg->Add(fullGoldGraph, "P"); 
  mg->Add(silverGraph, "P");
  mg->Add(goldGraph, "");  
  mg->Draw("AP");
  //mg->GetYaxis()->SetLimits(-1,3);
  //mg->GetXaxis()->SetLimits(0,80);
  mg->SetTitle("Gold and Silver Peak Voltages Over Total Illumination Time");
  mg->GetXaxis()->SetTitle("Time Illuminated (days)"); 
  mg->GetYaxis()->SetTitle("Peak Voltage (Volts)");
  auto legend = new TLegend(0.57,0.65,0.9,0.9);
  legend->SetNColumns(1);
  legend->SetTextSize(0.025);
  //legend->AddEntry(silverGraph->GetFunction("silver_inv"), "Fit p_{0} + #frac{p_{1}}{x+p_{2}}", "L");
  //legend->AddEntry(silverGraph, "Silver Peak Voltage", "P");
  //legend->AddEntry(silverGraph->GetFunction("silver_combo"), "Fit p_{0} + #frac{p_{1}}{p_{2}+x} + #frac{p_{3}}{x^{2}+p_{4}x+p_{5}}", "L");
  //legend->AddEntry(fullGoldGraph, "Gold Peak Voltage", "P");  
  //legend->AddEntry(silverGraph->GetFunction("silver_exp"), "Fit p_{0} + e^{p_{1}+p_{2}x}", "L");
  legend->AddEntry(fullGoldGraph->GetFunction("full_gold_combo"), "Fit including all points");
  legend->AddEntry(goldGraph->GetFunction("gold_combo"), "Fit excluding points after 5 days");
  legend->AddEntry(silverGraph, "Silver Peak Voltage", "P");
  legend->AddEntry(fullGoldGraph, "Gold Peak Voltage", "P");  
  legend->Draw();  
  
}