#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TTree.h"
#include "TClonesArray.h"

// ROOT includes
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TText.h"
#include "Riostream.h" //????????

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <TStyle.h>

// pt&y matrix
//--------------------------------
//  ptmax | h14 | h24 | h34 | h44 |
//     pt | h13 | h23 | h33 | h34 |
//     pt | h12 | h22 | h32 | h42 |
//  ptmin | h11 | h21 | h31 | h41 |
//--------------------------------
//          ymin   y     y    ymax

// pt&y matrix rebin
//--------------------------------------------------------

// ptmax | h18 | h28 | h38 | h48 | h58 | h68 | h78 | h88 |
//    pt | h17 | h27 | h37 | h47 | h57 | h67 | h77 | h87 |
//    pt | h16 | h26 | h36 | h46 | h56 | h66 | h76 | h86 |
//    pt | h15 | h25 | h35 | h45 | h55 | h65 | h75 | h85 |
//    pt | h14 | h24 | h34 | h44 | h54 | h64 | h74 | h84 |
//    pt | h13 | h23 | h33 | h34 | h53 | h63 | h73 | h83 |
//    pt | h12 | h22 | h32 | h42 | h52 | h62 | h72 | h82 |
// ptmin | h11 | h21 | h31 | h41 | h51 | h61 | h71 | h81 |
//--------------------------------------------------------
//         ymin   y     y     y     y     y     y    ymax

using namespace std;

// ########################################

   // таблица эффективностей 1/wi для реакции С + Сu(4ГэВ)

    Float_t coeff_cc45[8][8] =
       {{43.87014563, 44.68564356, 29.69687815, 37.50973236, 32.09756098, 38.54035568, 48.93617021, 68.71345029},
        {45.1439523, 40.45665172, 28.46015424, 31.29225352, 27.34043887, 28.98667602, 26.54188144, 47.16516517},
        {45.95238095, 31.89371981, 25.09760766, 21.8879273, 26.09016393, 21.46699752, 22.72571429, 25.16275168},
        {53.45302548, 33.18058252, 26.76818623, 23.98576676, 24.09373458, 20.46569855, 19.67382199, 21.5203252},
        {103.4935567, 45.76467269, 29.26318245, 23.94635545, 23.38011696, 20.66933467, 17.72552977, 21.13488844},
        {253.9963899, 83.77008653, 40.03873518, 26.47504404, 27.92797574, 22.80413223, 17.44200627, 32.13381995},
        {1667.974359, 161.3671875, 64.2936747, 41.15876777, 31.0304679, 27.84992785, 20.30896552, 42.57281553},
        {8399.4, 369.3235294, 133.3386243, 42.04068522, 33.89705882, 39.23293173, 22.84012539, 246.1875}};

   Float_t coeff_ccu45[8][8] =
       {{97.12665198, 56.82994652, 57.86950549, 35.4095064, 57.33065811, 35.16511628, 103.4481605, 71.16149068},
        {54.57709924, 43.81280521, 50.86515913, 41.20864947, 30.69267998, 30.58196135, 33.74306688, 52.48613678},
        {80.94619799, 44.71141553, 39.10123043, 34.12915327, 35.2310559, 24.71479959, 27.19510703, 64.09706546},
        {91.95371669, 46.83424908, 33.99166667, 30.61150401, 28.70030426, 29.79841897, 18.23875753, 25.86519337},
        {170.7622699, 83.02379286, 48.84844939, 33.00271862, 31.94380252, 27.73352941, 28.78217822, 24.13535589},
        {538.6542553, 110.4280936, 59.91189427, 46.05302403, 36.78336079, 25.91329011, 26.26190476, 27.41150442},
        {1096.711864, 217.8387909, 96.67179487, 54.04830918, 42.62165605, 25.89449541, 34.01443299, 63.9453125},
        {3081.238095, 755.8714286, 180.1823204, 74.74251497, 45.28301887, 29.7232376, 35.38461538, 69.76363636}};

   Float_t coeff_cal45[8][8] =
       {{62.01152369, 42.59747292, 49.45921053, 32.67639015, 39.00907029, 43.95903955, 29.71157324, 76.23291925},
        {48.57708049, 34.41529526, 29.07987048, 22.52075142, 34.45781466, 30.10584958, 29.48414795, 26.5572264},
        {51.66528681, 33.9184, 31.03726708, 26.29956332, 25.848, 20.85287659, 20.30206957, 22.26435935},
        {70.21769231, 39.63477461, 32.41064453, 24.23783573, 23.39830508, 24.782678, 20.18359375, 21.2375317},
        {129.6630952, 58.0438942, 34.84841418, 27.81255061, 25.57100466, 24.18304732, 2.665061898, 18.10318792},
        {295.3516129, 95.40528634, 55.30353982, 33.0489426, 32.2611997, 29.77981651, 19.67948718, 29.29912664},
        {1101.461538, 201.832021, 83.48019017, 45.8488121, 31.35707121, 22.31358529, 25.11229135, 46.69430052},
        {3974.928571, 795.4915254, 165.3459459, 78.57333333, 40.65411765, 28.67901235, 24.02824859, 57.36111111}};

   Float_t coeff_cpb45[8][8] =
       {{89.19794721, 51.46024322, 70.86481802, 48.64728682, 47.24099723, 34.25813692, 45.31516937, 35.94761905},
        {69.99224205, 52.77220077, 39.49669749, 47.15996578, 45.89859155, 52.99758454, 37.52858399, 34.60743322},
        {90.57813834, 60.93068683, 46.08560053, 43.03746594, 44.44110672, 37.01104972, 37.31753948, 39.63498623},
        {104.6717268, 63.66730159, 45.81899299, 39.97699005, 43.5177026, 31.82094834, 29.51099707, 49.62526767},
        {223.612069, 99.0, 54.18856767, 45.65893417, 37.18730554, 34.8575804, 29.14485981, 23.37485971},
        {554.8374384, 175.1385435, 94.48266297, 54.68380213, 43.71598808, 45.72172619, 24.18318318, 30.7675},
        {2390.047619, 295.4744027, 130.3418605, 61.66666667, 51.17391304, 50.69634703, 22.73735955, 67.31666667},
        {5797.909091, 1152.217391, 151.9722222, 99.67460317, 51.78678679, 33.07668712, 35.6, 147.1481481}};

   // ########################################

   //// end
   // pt & y intervals
   Float_t pt_intervals[4][2] = {{0.1, 0.3}, {0.3, 0.5}, {0.5, 0.75}, {0.75, 1.05}};   // 4;
   Float_t y_intervals[4][2] = {{1.2, 1.45}, {1.45, 1.65}, {1.65, 1.85}, {1.85, 2.1}}; // 4

   // ################ rebin by pt and rapidity(y)
   // 4.0 from MK
   Float_t pt_intervals_rebin[8][2] = {{0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.625}, {0.625, 0.75}, {0.75, 0.9}, {0.9, 1.05}}; // 4

   Float_t y_intervals_rebin[8][2] = {{1.2, 1.325}, {1.325, 1.45}, {1.45, 1.55}, {1.55, 1.65}, {1.65, 1.75}, {1.75, 1.85}, {1.85, 1.975}, {1.975, 2.1}}; // 4

   Float_t pt_intervals_hist40_rebin[9] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.625, 0.75, 0.9, 1.05};     // 4
   Float_t y_intervals_hist40_rebin[9] = {1.2, 1.325, 1.45, 1.55, 1.65, 1.75, 1.85, 1.975, 2.1}; // 4
   Float_t pt_intervals_hist[5] = {0.1, 0.3, 0.5, 0.75, 1.05};                                   // 4
   Float_t y_intervals_hist[5] = {1.2, 1.45, 1.65, 1.85, 2.1};              
   //################# rebin
   Float_t pt_intervals45_rebin[8][2]={{0.1,0.2},  {0.2,0.3},   {0.3,0.4}, {0.4,0.5},
				       {0.5,0.625},{0.625,0.75},{0.75,0.9},{0.9,1.05}};//45 NOTE
   
   Float_t y_intervals45_rebin[8][2]= {{1.25,1.375},{1.375,1.50}, {1.50,1.60},  {1.60,1.70},
				       {1.70,1.80}, {1.80, 1.90}, {1.90,2.025}, {2.025,2.15}};//45

   Float_t pt_intervals_hist45_rebin[9]={0.1, 0.2, 0.3, 0.4, 0.5, 0.625, 0.75, 0.9, 1.05}; //45

   Float_t y_intervals_hist45_rebin[9]={1.25, 1.375, 1.5, 1.6, 1.7, 1.8, 1.9, 2.025, 2.15}; //45                        // 4
   Float_t width_bins_pt[9];
   Float_t width_bins_y[9];
   Int_t Nbins = 36;// число разбиений для массовых гистограмм
   Double_t min_mass=1.078;
   Double_t max_mass=1.18;
   // ################ end rebin
   // ################
   Int_t bins = 4;
   Int_t bins_rebin = 8;

  void mass_L045(TString infile = "out_4.5GeV_CPb_pi3_p4_nMC2000_rebin.root"){

 // ##Read Tree - start
   TFile *file = new TFile(infile, "READ");
   TTree *tree = (TTree *)file->Get("data_inf");
   TH2D *histFit_p = (TH2D *)file->Get("hEff1dow2d_rebin");

   Double_t pts, ys, ms;

   tree->SetBranchAddress("data_mass", &ms);
   tree->SetBranchAddress("data_pt", &pts);
   tree->SetBranchAddress("data_y", &ys);

  Int_t pt_i, y_i;
   // коэфефициенты
   Float_t w_cc, w_cu, w_cal, w_cpb;
   TObjArray *Hlist = new TObjArray();

   // делаем 4 массовые гистограммы для областей pt и 4 массовые для областей y;
   TH1F *sig_mass_y_data_weight[4];
   TH1F *sig_mass_pt_data_weight[4];

   TString htitle = "";
   TString hname = "";

   // ######## for rebin делаем 64 гистаграммы для данных
   TH1F *sig_data_rebin[8][8];
   TH1F *sig_data_rebin_clon[8][8];
   TH1F *sig_data_rebin_clon_extr[8][8];

   // ###################### rebin
   for (Int_t i = 0; i < bins_rebin; i++)
   {
      for (Int_t j = 0; j < bins_rebin; j++)
      {
         hname.Form("hdata_rebin%d%d", i + 1, j + 1);
         htitle.Form("%3.2f<y<%3.2f, %3.2f<pt<%3.2f", y_intervals45_rebin[i][0], y_intervals45_rebin[i][1], pt_intervals45_rebin[j][0], pt_intervals45_rebin[j][1]);
         sig_data_rebin[i][j] = new TH1F(hname, "",  Nbins, min_mass, max_mass);
         sig_data_rebin[i][j]->SetTitle(htitle);
         sig_data_rebin[i][j]->SetMarkerStyle(20);
         Hlist->Add(sig_data_rebin[i][j]);
      }
   }

 // wцикл по событиям - начало
   Int_t events = tree->GetEntries();
   for (Int_t ev = 0; ev < events; ev++)
   {
      tree->GetEntry(ev);

      if ((ev % 10000) == 0)
         cout << " Event= " << ev << "/" << events << endl;
      if (pts < 0.1 && pts > 1.05 || ys < 1.2 && ys > 2.1)
       continue;

    Int_t i_sig = -999;
    Int_t j_sig = -999;
 // if(pts < 0.1 && pts > 1.05 || ys < 1.25 && ys > 2.15) continue;
   //определяем в какую ячейку попадает событие в данных по матрице 8x8 по pt & y          
  if(ys>=1.25 && ys<=1.38){
    y_i=0;
   }else if(ys>=1.38 && ys<=1.50){
    y_i=1;
   }else if(ys>=1.50 && ys<=1.60){
    y_i=2;
   }else if(ys>=1.60 && ys<=1.70){
    y_i=3;
   }else if(ys>=1.70 && ys<=1.80){
    y_i=4;
   }else if(ys>=1.80 && ys<=1.90){
    y_i=5;
   }else if(ys>=1.90 && ys<=2.03){
    y_i=6;
   }else if(ys>=2.03 && ys<=2.15){
    y_i=7;
   }else{
    continue;
   }
    
   if(pts>=0.1 && pts<=0.2){
    pt_i=0;
   }else if(pts>=0.2 && pts<=0.3){
    pt_i=1;
   }else if(pts>=0.3 && pts<=0.4){
    pt_i=2;
   }else if(pts>=0.4 && pts<=0.5){
    pt_i=3;
   }else if(pts>=0.5 && pts<=0.62){
    pt_i=4;
   }else if(pts>=0.62 && pts<=0.75){
    pt_i=5;
   }else if(pts>=0.75 && pts<=0.90){
    pt_i=6;
   }else if(pts>=0.90 && pts<=1.05){
    pt_i=7;
   }else{
   continue;
    }
     // w_cc = coeff_cc45[pt_i][y_i];
     // w_cal = coeff_cal45[pt_i][y_i];
      w_cpb = coeff_cpb45[pt_i][y_i];


    //  if (w_cc > 100)continue;
   // if (w_cal > 100)continue;
   if (w_cpb > 100)continue;
         // fill mass spectras in pt&y intervals
         for (Int_t pti = 0; pti < 8; pti++)
         {
            if (pts >= pt_intervals45_rebin[pti][0] && pts <= pt_intervals45_rebin[pti][1])
            {
               i_sig = pti + 1;
               break;
            }
         }
         for (Int_t yj = 0; yj < 8; yj++)
         {
            if (ys >= y_intervals45_rebin[yj][0] && ys <= y_intervals45_rebin[yj][1])
            {
               j_sig = yj + 1;

               break;
            }
         }
         if (i_sig > 0 && j_sig > 0)
         {
            //sig_data_rebin[j_sig - 1][i_sig - 1]->Fill(ms, w_cal);
            sig_data_rebin[j_sig - 1][i_sig - 1]->Fill(ms, w_cpb);
         }
   }

Float_t Val[8][8];
Float_t Val_Sum;

for (Int_t k = 0; k < 8; k++)
{
         for (Int_t kk = 0; kk < 8; kk++)
         {
            width_bins_pt[k] = pt_intervals45_rebin[k][1] - pt_intervals45_rebin[k][0];
            width_bins_y[kk] = y_intervals45_rebin[kk][1] - y_intervals45_rebin[kk][0];

            Val[k][kk] = width_bins_pt[k] * width_bins_y[kk];
            Val_Sum += Val[k][kk];
            cout << "  " << width_bins_pt[k] << " ,  " << width_bins_y[kk] << " , " << Val[k][kk] <<endl;
         }
}
//cout << " sum = " <<  Val_Sum << endl;
// TObjArray *Hlist = new TObjArray();
Double_t s, err, ds, derr;
Double_t ent;
for (Int_t i = 0; i < bins_rebin; i++)
{
         for (Int_t j = 0; j < bins_rebin; j++)
         {
            sig_data_rebin_clon[i][j] = (TH1F *)sig_data_rebin[i][j]->Clone("sig_data_rebin_clon");
           // sig_data_rebin_clon_CAl45[i][j] = (TH1F *)sig_data_rebin[i][j]->Clone("sig_data_rebin_clon_CAl45");
            //Hlist->Add(sig_data_rebin_clon[i][j]);

            for (Int_t k = 0; k < sig_data_rebin[i][j]->GetNbinsX(); k++)
            {

               ent = sig_data_rebin[i][j]->GetEntries();
               if (ent < 6)
                  continue;

               s = sig_data_rebin[i][j]->GetBinContent(k + 1);
               err = sig_data_rebin[i][j]->GetBinError(k + 1);
               if (s == 0 && err == 0)
                  continue;
               // cout << "  " << s << " , " << err << endl;
               //   ds = s / Val[i][j];
               //  derr = err / Val[i][j];
               ds = s * Val_Sum/Val[i][j];
                derr = err * Val_Sum/Val[i][j];
               //ds = s * Val[i][j] / Val_Sum;
              // derr = err * Val[i][j] / Val_Sum;

               if (ds == 0 && derr == 0)
                  continue;
               // cout << "  " << ds << " , " << derr << endl;
               // if(i>=j){
               sig_data_rebin_clon[i][j]->SetBinContent(k + 1, ds);
               sig_data_rebin_clon[i][j]->SetBinError(k + 1, derr);
               //}
           } 
           
         }
}
 

 //###################Procedure extrapolation##################################
//wextrp_con[4];

for(Int_t ik=0;ik<bins_rebin; ik++){
      for(Int_t jk=0;jk<bins_rebin; jk++){

      sig_data_rebin_clon_extr[ik][jk]=(TH1F*) sig_data_rebin_clon[ik][jk]->Clone("sig_data_rebin_clon_extr"); 
      //sig_data_rebin_clon_extr[ik][jk]->Reset();
      //   sig_data_rebin_clon[i][j]= (TH1F*) sig_data_rebin[i][j]->Clone("sig_data_rebin_clon"); 
     // Hlist->Add(sig_data_rebin_clon_extr[ik][jk]);
      }
   }
//Double_t wextrp_con_CCu[4]={1.493700503 ,1.976310226 ,1.426400599 ,1.024563627};

//Double_t wextrp_con[4]={1. ,1. ,1. ,1.};
//Double_t wextrp_con[4] ={2.106650699, 1.526689315, 1.071798112, 1.02657736};//C+Cu(4.5)
//Double_t wextrp_con[4]={2.1665, 1.5476, 1.0747, 1.0260};//C+C(4.5)
//Double_t wextrp_con[3]={2.14665062 ,1.258503572 ,1.074119933};//C+Al(4.5)
Double_t wextrp_con[4]={3.009995333,1.527522397,1.226355938,1.027285963};//C+Pb(4.5)



Double_t s_mass_data[8][8];
Double_t err_mass_data[8][8];
Double_t biny;
Double_t bin_err;
for (Int_t ki = 0; ki < bins_rebin; ki++)
{
      for (Int_t kj = 0; kj < bins_rebin; kj++)
      {
      for (Int_t h = 0; h < sig_data_rebin_clon[ki][kj]->GetNbinsX(); h++)
      {
               s_mass_data[ki][kj] = sig_data_rebin[ki][kj]->GetBinContent(h + 1);
               err_mass_data[ki][kj] = sig_data_rebin[ki][kj]->GetBinError(h + 1);

               if(s_mass_data[ki][kj]==0 && err_mass_data[ki][kj]==0) continue;

               cout << " getcontent for extrapolate " << s_mass_data[ki][kj] << endl;
               cout << "Error for extrapolate " << err_mass_data[ki][kj] << endl;


               if (ki == 0)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
               }
               else if (ki == 1)
               {
                  biny = wextrp_con[1] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[1] * err_mass_data[ki][kj];
               }
               else if (ki == 2)
               {
                  biny = wextrp_con[2] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[2] * err_mass_data[ki][kj];
               } 
              else if (ki == 7)
               {
                  biny = wextrp_con[3] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[3] * err_mass_data[ki][kj];
               }
               else
               {
                  biny = s_mass_data[ki][kj];
                  bin_err = err_mass_data[ki][kj];
               }

               sig_data_rebin_clon_extr[ki][kj]->SetBinContent(h + 1, biny);
               sig_data_rebin_clon_extr[ki][kj]->SetBinError(h + 1, bin_err);
      }
      }
}
   TH1F *summ_hist_y1 = new TH1F("summ_hist_y1", "1.20 < y < 1.45", Nbins, min_mass, max_mass);
   summ_hist_y1->SetMarkerStyle(20);
   TH1F *summ_hist_y2 = new TH1F("summ_hist_y2", " ", Nbins, min_mass, max_mass);
   summ_hist_y2->SetMarkerStyle(20);
   TH1F *summ_hist_y3 = new TH1F("summ_hist_y3", "1.45 < y < 1.65", Nbins, min_mass, max_mass);
   summ_hist_y3->SetMarkerStyle(20);
   TH1F *summ_hist_y4 = new TH1F("summ_hist_y4", "", Nbins, min_mass, max_mass);
   summ_hist_y4->SetMarkerStyle(20);
   TH1F *summ_hist_y5 = new TH1F("summ_hist_y5", "1.65 < y < 1.85", Nbins, min_mass, max_mass);
   summ_hist_y5->SetMarkerStyle(20);
   TH1F *summ_hist_y6 = new TH1F("summ_hist_y6", "", Nbins, min_mass, max_mass);
   summ_hist_y6->SetMarkerStyle(20);
   TH1F *summ_hist_y7 = new TH1F("summ_hist_y7", "1.85 < y < 2.1", Nbins, min_mass, max_mass);
   summ_hist_y7->SetMarkerStyle(20);
   TH1F *summ_hist_y8 = new TH1F("summ_hist_y8", "", Nbins, min_mass, max_mass);
   summ_hist_y8->SetMarkerStyle(20);

   //===================== Суммируем строки =============================================

   TH1F *summ_hist_p1 = new TH1F("summ_hist_p1", "0.10 < pt < 0.30", Nbins, min_mass, max_mass);
   summ_hist_p1->SetMarkerStyle(20);
   TH1F *summ_hist_p2 = new TH1F("summ_hist_p2", "", Nbins, min_mass, max_mass);
   summ_hist_p2->SetMarkerStyle(20);
   TH1F *summ_hist_p3 = new TH1F("summ_hist_p3", "0.30 < pt < 0.50", Nbins, min_mass, max_mass);
   summ_hist_p3->SetMarkerStyle(20);
   TH1F *summ_hist_p4 = new TH1F("summ_hist_p4", "", Nbins, min_mass, max_mass);
   summ_hist_p4->SetMarkerStyle(20);
   TH1F *summ_hist_p5 = new TH1F("summ_hist_p5", "0.50 < pt < 0.75", Nbins, min_mass, max_mass);
   summ_hist_p5->SetMarkerStyle(20);
   TH1F *summ_hist_p6 = new TH1F("summ_hist_p6", "", Nbins, min_mass, max_mass);
   summ_hist_p6->SetMarkerStyle(20);
   TH1F *summ_hist_p7 = new TH1F("summ_hist_p7", "0.75 < pt < 1.05", Nbins, min_mass, max_mass);
   summ_hist_p7->SetMarkerStyle(20);
   TH1F *summ_hist_p8 = new TH1F("summ_hist_p8", "", Nbins, min_mass, max_mass);
   summ_hist_p8->SetMarkerStyle(20);

   for (Int_t j0 = 0; j0 < 8; j0++)
   {
      summ_hist_y1->Add(sig_data_rebin_clon_extr[0][j0], 1);
   }
   //}
   // for (Int_t i1=1; i1<2; i1++){
   for (Int_t j1 = 0; j1 < 8; j1++)
   {
      summ_hist_y2->Add(sig_data_rebin_clon_extr[1][j1], 1);
   }
   //}
   // for (Int_t i2=2; i2<3; i2++){
   for (Int_t j2 = 0; j2 < 8; j2++)
   {
      summ_hist_y3->Add(sig_data_rebin_clon_extr[2][j2], 1);
   }
   //}
   // for (Int_t i3=3; i3<4; i3++){
   for (Int_t j3 = 0; j3 < 8; j3++)
   {
      summ_hist_y4->Add(sig_data_rebin_clon_extr[3][j3], 1);
   }
   //}
   // for (Int_t i4=4; i4<5; i4++){
   for (Int_t j4 = 0; j4 < 8; j4++)
   {
      summ_hist_y5->Add(sig_data_rebin_clon_extr[4][j4], 1);
   }
   //}
   // for (Int_t i5=5; i5<6; i5++){
   for (Int_t j5 = 0; j5 < 8; j5++)
   {
      summ_hist_y6->Add(sig_data_rebin_clon_extr[5][j5], 1);
   }
   //}
   // for (Int_t i6=6; i6<7; i6++){
   for (Int_t j6 = 0; j6 < 8; j6++)
   {
      summ_hist_y7->Add(sig_data_rebin_clon_extr[6][j6], 1);
   }
   //}
   // for (Int_t i7=7; i7<8; i7++){
   for (Int_t j7 = 0; j7 < 8; j7++)
   {
      summ_hist_y8->Add(sig_data_rebin_clon_extr[7][j7], 1);
   }
   //}
     for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      // for (Int_t jj0=0; jj0<1; jj0++){
      summ_hist_p1->Add(sig_data_rebin_clon_extr[ii0][0], 1);
      //  }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      // for (Int_t jj0=1; jj0<2; jj0++){
      summ_hist_p2->Add(sig_data_rebin_clon_extr[ii0][1], 1);
      //  }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      // for (Int_t jj0=2; jj0<3; jj0++){
      summ_hist_p3->Add(sig_data_rebin_clon_extr[ii0][2], 1);
      //  }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      // for (Int_t jj0=3; jj0<4; jj0++){
      summ_hist_p4->Add(sig_data_rebin_clon_extr[ii0][3], 1);
      // }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      // for (Int_t jj0=4; jj0<5; jj0++){
      summ_hist_p5->Add(sig_data_rebin_clon_extr[ii0][4], 1);
      // }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      // for (Int_t jj0=5; jj0<6; jj0++){
      summ_hist_p6->Add(sig_data_rebin_clon_extr[ii0][5], 1);
      // }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      //  for (Int_t jj0=6; jj0<7; jj0++){
      summ_hist_p7->Add(sig_data_rebin_clon_extr[ii0][6], 1);
      //  }
   }
   for (Int_t ii0 = 0; ii0 < 8; ii0++)
   {
      //  for (Int_t jj0=7; jj0<7; jj0++){
      summ_hist_p8->Add(sig_data_rebin_clon_extr[ii0][7], 1);
      // }
   }

    ////Cуммировать парами 8 гистограмм по быстроте и получить 4 гистограммы по быстроте

   TH1F *summ_4y_b1 = (TH1F *)summ_hist_y1->Clone("summ_4y_b1");
   summ_4y_b1->Add(summ_hist_y2, 1);
   TH1F *summ_4y_b2 = (TH1F *)summ_hist_y3->Clone("summ_4y_b2");
   summ_4y_b2->Add(summ_hist_y4, 1);
   TH1F *summ_4y_b3 = (TH1F *)summ_hist_y5->Clone("summ_4y_b3");
   summ_4y_b3->Add(summ_hist_y6, 1);
   TH1F *summ_4y_b4 = (TH1F *)summ_hist_y7->Clone("summ_4y_b4");
   summ_4y_b4->Add(summ_hist_y8, 1);
   //=================Sum Hist суммируем строчки//
   Hlist->Add(summ_4y_b1);
   Hlist->Add(summ_4y_b2);
   Hlist->Add(summ_4y_b3);
   Hlist->Add(summ_4y_b4);

   ////Cуммировать парами 8 гистограмм по р_Т и получить 4 гистограммы по р_Т.
   TH1F *summ_4pt_b1 = (TH1F *)summ_hist_p1->Clone("summ_4pt_b1");
   summ_4pt_b1->Add(summ_hist_p2, 1);
   TH1F *summ_4pt_b2 = (TH1F *)summ_hist_p3->Clone("summ_4pt_b2");
   summ_4pt_b2->Add(summ_hist_p4, 1);
   TH1F *summ_4pt_b3 = (TH1F *)summ_hist_p5->Clone("summ_4pt_b3");
   summ_4pt_b3->Add(summ_hist_p6, 1);
   TH1F *summ_4pt_b4 = (TH1F *)summ_hist_p7->Clone("summ_4pt_b4");
   summ_4pt_b4->Add(summ_hist_p8, 1);

   Hlist->Add(summ_4pt_b1);
   Hlist->Add(summ_4pt_b2);
   Hlist->Add(summ_4pt_b3);
   Hlist->Add(summ_4pt_b4); 


   TH1F *widht_4pt_b1 = (TH1F *)summ_4pt_b1->Clone("widht_4pt_b1");
   TH1F *widht_4pt_b2 = (TH1F *)summ_4pt_b2->Clone("widht_4pt_b2");
   TH1F *widht_4pt_b3 = (TH1F *)summ_4pt_b3->Clone("widht_4pt_b3");
   TH1F *widht_4pt_b4 = (TH1F *)summ_4pt_b4->Clone("widht_4pt_b4");
   
   TH1F *hist_sum_by_pt_1 = (TH1F *)widht_4pt_b1->Clone("hist_sum_by_pt_1");
   hist_sum_by_pt_1->Add(widht_4pt_b2, 1);
   TH1F *hist_sum_by_pt_3 = (TH1F *)hist_sum_by_pt_1->Clone("hist_sum_by_pt_3");
   hist_sum_by_pt_3->Add(widht_4pt_b3, 1);
   TH1F *hist_sum_by_pt = (TH1F *)hist_sum_by_pt_3->Clone("hist_sum_by_pt");
   hist_sum_by_pt->Add(widht_4pt_b4, 1);

   Hlist->Add(hist_sum_by_pt);
   //=============================от 4 к 1 гистограмме по быстроте=====//
   TH1F *widht_4y_b1 = (TH1F *)summ_4y_b1->Clone("widht_4y_b1");
   TH1F *widht_4y_b2 = (TH1F *)summ_4y_b2->Clone("widht_4y_b3");
   TH1F *widht_4y_b3 = (TH1F *)summ_4y_b3->Clone("widht_4y_b2"); 
   TH1F *widht_4y_b4 = (TH1F *)summ_4y_b4->Clone("widht_4y_b4");


   TH1F *hist_sum_by_y_1 = (TH1F *)widht_4y_b1->Clone("hist_sum_by_y_1");
   hist_sum_by_y_1->Add(widht_4y_b2, 1);
   TH1F *hist_sum_by_y_3 = (TH1F *)hist_sum_by_y_1->Clone("hist_sum_by_y_3");
   hist_sum_by_y_3->Add(widht_4y_b3, 1);
   TH1F *hist_sum_by_y = (TH1F *)hist_sum_by_y_3->Clone("hist_sum_by_y");
   hist_sum_by_y->Add(widht_4y_b4, 1);

   Hlist->Add(hist_sum_by_y);
     
   TFile *out = new TFile("outFile_test_29082023_CPb45Gev_p2.root", "RECREATE");
   Hlist->Write();
   out->Write();
   cout << " WRITE FILE  " << endl;
   out->Close();
   }