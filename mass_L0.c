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
Float_t coeff_ccu[8][8] =
    {{0.0201, 0.0237, 0.0222, 0.0289, 0.0398, 0.0258, 0.0570, 0.0502},
     {0.0188, 0.0300, 0.0342, 0.0438, 0.0355, 0.0443, 0.0440, 0.0483},
     {0.0180, 0.0197, 0.0289, 0.0436, 0.0336, 0.0534, 0.0386, 0.0428},
     {0.0171, 0.0156, 0.0222, 0.0392, 0.0468, 0.0427, 0.0627, 0.0419},
     {0.0111, 0.0134, 0.0126, 0.0335, 0.0296, 0.0449, 0.0419, 0.0367},
     {0.0061, 0.0086, 0.0070, 0.0219, 0.0300, 0.0379, 0.0454, 0.0090},
     {0.0020, 0.0013, 0.0018, 0.0126, 0.0353, 0.0277, 0.0153, 0.0184},
     {0.0002, 0.0002, 0.0004, 0.0135, 0.0138, 0.0146, 0.0257, 0.0235}};

Float_t coeff_cc[8][8] =
    {{0.0308, 0.0362, 0.0329, 0.0486, 0.0488, 0.0520, 0.0667, 0.0618},
     {0.0255, 0.0329, 0.0375, 0.0482, 0.0516, 0.0378, 0.0481, 0.0572},
     {0.0236, 0.0355, 0.0442, 0.0582, 0.0553, 0.0530, 0.0516, 0.0334},
     {0.0141, 0.0290, 0.0420, 0.0520, 0.0580, 0.0582, 0.0421, 0.0546},
     {0.0058, 0.0172, 0.0313, 0.0424, 0.0478, 0.0518, 0.0487, 0.0291},
     {0.0020, 0.0082, 0.0216, 0.0304, 0.0414, 0.0369, 0.0505, 0.0251},
     {0.0001, 0.0031, 0.0093, 0.0185, 0.0264, 0.0449, 0.0335, 0.0189},
     {0.0002, 0.0007, 0.0040, 0.0150, 0.0225, 0.0231, 0.0248, 0.0083}};

Float_t coeff_cal[8][8] = {{0.0289, 0.0243, 0.0293, 0.0422, 0.0392, 0.0438, 0.0348, 0.0241},
                           {0.0204, 0.0344, 0.0340, 0.0349, 0.0393, 0.0495, 0.0491, 0.0708},
                           {0.0193, 0.0280, 0.0353, 0.0471, 0.0535, 0.0555, 0.0452, 0.0445},
                           {0.0104, 0.0208, 0.0361, 0.0453, 0.0484, 0.0500, 0.0548, 0.0297},
                           {0.0051, 0.0163, 0.0238, 0.0306, 0.0426, 0.0435, 0.0440, 0.0332},
                           {0.0012, 0.0072, 0.0137, 0.0298, 0.0345, 0.0439, 0.0429, 0.0258},
                           {0.0004, 0.0025, 0.0113, 0.0123, 0.0235, 0.0422, 0.0183, 0.0105},
                           {0.0002, 0.0002, 0.0050, 0.0077, 0.0180, 0.0285, 0.0100, 0.0057}};

Float_t coeff_cpb[8][8] =
    {{0.016632692, 0.025368145, 0.027722825, 0.032728373, 0.037890596, 0.040889558, 0.049039459, 0.044258648},
     {0.012684313, 0.017328327, 0.025700409, 0.031340605, 0.037079089, 0.028043605, 0.048419574, 0.019942316},
     {0.007810197, 0.014848752, 0.023728113, 0.028759786, 0.032839754, 0.042692578, 0.023418559, 0.040286482},
     {0.00345458, 0.008853897, 0.015010966, 0.029347635, 0.028606757, 0.039418994, 0.022441539, 0.035868522},
     {0.001486216, 0.005259155, 0.013527296, 0.020209505, 0.028857753, 0.024543695, 0.041302098, 0.020224719},
     {0.000402564, 0.00163912, 0.005574936, 0.011643627, 0.017726734, 0.019976383, 0.029585799, 0.005619639},
     {0.000480392, 0.000918153, 0.002997708, 0.006536482, 0.019057172, 0.030728376, 0.007233796, 0.007472826}};

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
Float_t y_intervals_hist[5] = {1.2, 1.45, 1.65, 1.85, 2.1};                                   // 4
Float_t width_bins_pt[9];
Float_t width_bins_y[9];
// ################ end rebin
// ################
Int_t bins = 4;
Int_t bins_rebin = 8;
Int_t Nbins = 39; // число разбиений для массовых гистограмм
Double_t min_mass = 1.078;
Double_t max_mass = 1.18;
   TString htitle = "";
   TString hname = "";
   TH1F *sig_data_rebin[8][8];
   TH1F *sig_data_rebin_clon[8][8];
   TH1F *sig_data_rebin_clon_extr[8][8];
   Double_t pts, ys, ms;
   Int_t pt_i, y_i;
   // коэфефициенты
   Float_t w_cc, w_cu, w_cal, w_cpb;

 // void mass_L0_CC(TString infile = "out_4.0GeV_CC_pi3_p4_nMC2000_rebin.root")
    //void mass_L0_CC(TString infile = "out_4.0GeV_CC_pi3_p4_nMC2000_rebin.root")
   void mass_L0_CC(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CC_pi3_p4_nMC2000_periodII_rebin.root")
   // void mass_L0_CC(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CC_pi3_p4_nMC2000_periodII_rebin.root")


   {
      // ##Read Tree - start
      TFile *file = new TFile(infile, "READ");
      TTree *tree = (TTree *)file->Get("data_inf");     
    
      tree->SetBranchAddress("data_mass", &ms);
      tree->SetBranchAddress("data_pt", &pts);
      tree->SetBranchAddress("data_y", &ys);
    
      TObjArray *Hlist = new TObjArray();
      // ###################### rebin
      for (Int_t i = 0; i < bins_rebin; i++)
      {
         for (Int_t j = 0; j < bins_rebin; j++)
         {
            hname.Form("hdata_rebin%d%d", i + 1, j + 1);
            htitle.Form("%3.2f<y<%3.2f, %3.2f<pt<%3.2f", y_intervals_rebin[i][0], y_intervals_rebin[i][1], pt_intervals_rebin[j][0], pt_intervals_rebin[j][1]);
            sig_data_rebin[i][j] = new TH1F(hname, "", Nbins, min_mass, max_mass);
            sig_data_rebin[i][j]->SetTitle(htitle);
            sig_data_rebin[i][j]->SetMarkerStyle(20);
            sig_data_rebin[i][j]->Sumw2();
           // Hlist->Add(sig_data_rebin[i][j]);
         }
      }

      // wцикл по событиям - начало
      Int_t events = tree->GetEntries();
      for (Int_t ev = 0; ev < events; ev++)
      {
         tree->GetEntry(ev);

         if ((ev % 10000) == 0)
            cout << " Event= " << ev << "/" << events << endl;
         // if (pts < 0.1 && pts > 1.05 || ys < 1.2 && ys > 2.1)
         //  continue;

         Int_t i_sig = -999;
         Int_t j_sig = -999;
         if (ys >= 1.2 && ys <= 1.33)
         {
            y_i = 0;
         }
         else if (ys >= 1.33 && ys <= 1.45)
         {
            y_i = 1;
         }
         else if (ys >= 1.45 && ys <= 1.55)
         {
            y_i = 2;
         }
         else if (ys >= 1.55 && ys <= 1.65)
         {
            y_i = 3;
         }
         else if (ys >= 1.65 && ys <= 1.75)
         {
            y_i = 4;
         }
         else if (ys >= 1.75 && ys <= 1.85)
         {
            y_i = 5;
         }
         else if (ys >= 1.85 && ys <= 1.98)
         {
            y_i = 6;
         }
         else if (ys >= 1.98 && ys <= 2.1)
         {
            y_i = 7;
         }
         /*  else
           {
              continue;
           }*/

         if (pts >= 0.1 && pts <= 0.2)
         {
            pt_i = 0;
         }
         else if (pts >= 0.2 && pts <= 0.3)
         {
            pt_i = 1;
         }
         else if (pts >= 0.3 && pts <= 0.4)
         {
            pt_i = 2;
         }
         else if (pts >= 0.4 && pts <= 0.5)
         {
            pt_i = 3;
         }
         else if (pts >= 0.5 && pts <= 0.62)
         {
            pt_i = 4;
         }
         else if (pts >= 0.62 && pts <= 0.75)
         {
            pt_i = 5;
         }
         else if (pts >= 0.75 && pts <= 0.90)
         {
            pt_i = 6;
         }
         else if (pts >= 0.90 && pts <= 1.05)
         {
            pt_i = 7;
         }
         /* else
          {
             continue;
          }*/
         // w_cu = 1. / coeff_ccu[pt_i][y_i];
         w_cc = 1. / coeff_cc[pt_i][y_i];
          if (w_cc > 100)
            continue;
         //  fill mass spectras in pt&y intervals
         for (Int_t pti = 0; pti < 8; pti++)
         {
            if (pts >= pt_intervals_rebin[pti][0] && pts <= pt_intervals_rebin[pti][1])
            {
               i_sig = pti + 1;
               break;
            }
         }
         for (Int_t yj = 0; yj < 8; yj++)
         {
            if (ys >= y_intervals_rebin[yj][0] && ys <= y_intervals_rebin[yj][1])
            {
               j_sig = yj + 1;

               break;
            }
         }
         if (i_sig > 0 && j_sig > 0)
         {
            sig_data_rebin[j_sig - 1][i_sig - 1]->Fill(ms, w_cc);
          
         }
}

Float_t Val[8][8];
Float_t Val_Sum;

for (Int_t k = 0; k < 8; k++)
{
         for (Int_t kk = 0; kk < 8; kk++)
         {
            width_bins_pt[k] = pt_intervals_rebin[k][1] - pt_intervals_rebin[k][0];
            width_bins_y[kk] = y_intervals_rebin[kk][1] - y_intervals_rebin[kk][0];

            Val[k][kk] = width_bins_pt[k] * width_bins_y[kk];
            Val_Sum += Val[k][kk];
           // cout << "  " << width_bins_pt[k] << " ,  " << width_bins_y[kk] << " , " << Val[k][kk] <<endl;
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
            // Hlist->Add(sig_data_rebin_clon[i][j]);

            for (Int_t k = 0; k < sig_data_rebin[i][j]->GetNbinsX(); k++)
            {
               ent = sig_data_rebin[i][j]->GetEntries();

               if(ent < 2)continue;
               // cout << "ent =  " << ent << endl;
               s = sig_data_rebin[i][j]->GetBinContent(k + 1);
               err = sig_data_rebin[i][j]->GetBinError(k + 1);
               if (s == 0 && err == 0)
                  continue;
           //    cout << "  " << s << " , " << err << endl;

               ds = s * Val[i][j] / Val_Sum;
               derr = err * Val[i][j] / Val_Sum;

               if (ds == 0 && derr == 0)
                  continue;
               //cout << "  " << ds << " , " << derr << endl;
               // if(i>=j){
               sig_data_rebin_clon[i][j]->SetBinContent(k + 1, ds);
               sig_data_rebin_clon[i][j]->SetBinError(k + 1, derr);
               //  }
               //}
            }
         }
}
 //###################Procedure extrapolation##################################
for (Int_t ik = 0; ik < bins_rebin; ik++)
{
         for (Int_t jk = 0; jk < bins_rebin; jk++)
         {

            sig_data_rebin_clon_extr[ik][jk] = (TH1F *)sig_data_rebin_clon[ik][jk]->Clone("sig_data_rebin_clon_extr");
         }
}
  Double_t wextrp_con[3] = {2.03056, 1.091547, 1.008}; // reaction C+C
 /*  Double_t s_mass_data;
   Double_t err_mass_data;
   Double_t biny;
   Double_t bin_err;*/
   //for (Int_t ki = 0; ki < bins_rebin; ki++)
   //{
    //  for (Int_t kj = 0; kj < bins_rebin; kj++)
      //{
  /* TH1F *s4y_b1 = (TH1F *)summ_4y_b1->Clone("s4y_b1");
   for (Int_t h = 0; h < summ_4y_b1->GetNbinsX(); h++)
   {
      s_mass_data = summ_4y_b1->GetBinContent(h + 1);
      err_mass_data = summ_4y_b1->GetBinError(h + 1);

      // if (s_mass_data[ki][kj] == 0 && err_mass_data[ki][kj] == 0)
      //  continue;
      biny = wextrp_con[0] * s_mass_data;
      bin_err = wextrp_con[0] * err_mass_data;

      s4y_b1->SetBinContent(h + 1, biny);
      s4y_b1->SetBinError(h + 1, bin_err);
   }
   // }
  // Hlist->Add(s4y_b1);
   Double_t s_mass_data_b2;
   Double_t err_mass_data_b2;
   Double_t biny_b2;
   Double_t bin_err_b2;
   // }
  TH1F *s4y_b2 = (TH1F *)summ_4y_b2->Clone("s4y_b2");
   for (Int_t h = 0; h < summ_4y_b2->GetNbinsX(); h++)
   {
      s_mass_data_b2 = summ_4y_b2->GetBinContent(h + 1);
      err_mass_data_b2 = summ_4y_b2->GetBinError(h + 1);

      // if (s_mass_data[ki][kj] == 0 && err_mass_data[ki][kj] == 0)
      //  continue;
      biny_b2 = wextrp_con[1] * s_mass_data_b2;
      bin_err_b2 = wextrp_con[1] * err_mass_data_b2;

      s4y_b2->SetBinContent(h + 1, biny_b2);
      s4y_b2->SetBinError(h + 1, bin_err_b2);
   }*/
   // }
  // Hlist->Add(s4y_b2);
 // Double_t wextrp_con[3] = {2.03056, 1.091547, 1.008}; // reaction C+C
//Double_t wextrp_con[4]={2.7016,1.8512,1.1270,1.0519};//reaction C+C
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

               if (s_mass_data[ki][kj] == 0 && err_mass_data[ki][kj] == 0)
                  continue;

               if (ki == 0)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
               }
               else if (ki == 1)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
               }
               else if (ki == 2)
               {
                  biny = wextrp_con[1] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[1] * err_mass_data[ki][kj];
               }
               else if (ki == 3)
               {
                  biny = wextrp_con[1] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[1] * err_mass_data[ki][kj];
               }
               else if (ki == 6)
               {
                  biny = wextrp_con[2] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[2] * err_mass_data[ki][kj];
               }
               else if (ki == 7)
               {
                  biny = wextrp_con[2] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[2] * err_mass_data[ki][kj];
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
/*Double_t s_mass_data[8][8];
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
}*/

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

    // }
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
/*    TH1F *widht_4pt_b1 = (TH1F *)s4y_b1->Clone("widht_4pt_b1");
   TH1F *widht_4pt_b2 = (TH1F *)s4y_b2->Clone("widht_4pt_b2");
   TH1F *widht_4pt_b3 = (TH1F *)summ_4pt_b3->Clone("widht_4pt_b3");
   TH1F *widht_4pt_b4 = (TH1F *)summ_4pt_b4->Clone("widht_4pt_b4");*/

   TH1F *hist_sum_by_pt_1 = (TH1F *)widht_4pt_b1->Clone("hist_sum_by_pt_1");
   hist_sum_by_pt_1->Add(widht_4pt_b2, 1);
   TH1F *hist_sum_by_pt_3 = (TH1F *)hist_sum_by_pt_1->Clone("hist_sum_by_pt_3");
   hist_sum_by_pt_3->Add(widht_4pt_b3, 1);
   TH1F *hist_sum_by_pt = (TH1F *)hist_sum_by_pt_3->Clone("hist_sum_by_pt");
   hist_sum_by_pt->Add(widht_4pt_b4, 1);
   //hist_sum_by_pt->Sumw2();
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
   //hist_sum_by_y->Sumw2();
   Hlist->Add(hist_sum_by_y);
 
  ////Cуммировать парами 8 гистограмм по быстроте и получить 4 гистограммы по быстроте

   TH1F *hpt_b1 = (TH1F *)summ_4pt_b1->Clone("hpt_b1");
   TH1F *hpt_b2 = (TH1F *)summ_4pt_b2->Clone("hpt_b2");
   TH1F *hpt_b3 = (TH1F *)summ_4pt_b3->Clone("hpt_b3");
   TH1F *hpt_b4 = (TH1F *)summ_4pt_b4->Clone("hpt_b4");
   //=================Sum Hist суммируем строчки//
   Double_t s_pt,s_pt2,s_pt3,s_pt4;
   Double_t err_pt,err_pt2,err_pt3,err_pt4;
   Double_t bin1,bin2,bin3,bin4;
   Double_t bin_err1,bin_err2,bin_err3,bin_err4;

   for (Int_t h = 0; h < summ_4pt_b1->GetNbinsX(); h++)
   {
      s_pt = summ_4pt_b1->GetBinContent(h + 1);
      err_pt = summ_4pt_b1->GetBinError(h + 1);

     // bin1 =0.3687*s_pt;Iper
      //bin_err1=0.3687*err_pt;Iper
     //     bin1 =0.631*s_pt;
    //  bin_err1=0.631*err_pt;

      hpt_b1->SetBinContent(h + 1, bin1);
      hpt_b1->SetBinError(h + 1, bin_err1);
   }
     for (Int_t h1 = 0; h1 < summ_4pt_b2->GetNbinsX(); h1++)
   {
      s_pt2 = summ_4pt_b2->GetBinContent(h1 + 1);
      err_pt2 = summ_4pt_b2->GetBinError(h1 + 1);

      bin2 =0.386*s_pt2;
      bin_err2=0.386*err_pt2;

      hpt_b2->SetBinContent(h1 + 1, bin2);
      hpt_b2->SetBinError(h1 + 1, bin_err2);
   }

    for (Int_t h2 = 0; h2< summ_4pt_b3->GetNbinsX(); h2++)
   {
      s_pt3 = summ_4pt_b3->GetBinContent(h2 + 1);
      err_pt3 = summ_4pt_b3->GetBinError(h2 + 1);

      bin3 =0.631*s_pt3;
      bin_err3=0.631*err_pt3;

      hpt_b3->SetBinContent(h2 + 1, bin3);
      hpt_b3->SetBinError(h2 + 1, bin_err3);
   }
    for (Int_t h3 = 0; h3 < summ_4pt_b4->GetNbinsX(); h3++)
   {
      s_pt4 = summ_4pt_b4->GetBinContent(h3 + 1);
      err_pt4 = summ_4pt_b4->GetBinError(h3 + 1);

      bin4 =0.571*s_pt4;
      bin_err4=0.571*err_pt4;

      hpt_b4->SetBinContent(h3 + 1, bin4);
      hpt_b4->SetBinError(h3 + 1, bin_err4);
   }

   Hlist->Add(hpt_b1);
   Hlist->Add(hpt_b2);
   Hlist->Add(hpt_b3);
   Hlist->Add(hpt_b4);
   //======================================суммирование диференцильное значение==========

  // TFile *out = new TFile("outFile_full_8bins_part_CC_periodI.root", "RECREATE");
   //TFile *out = new TFile("outFile_full_8bins_part_CC_periodII.root", "RECREATE");
  TFile *out = new TFile("outFile_full_8bins_part_CC_full_refit_perII.root", "RECREATE");
   Hlist->Write();
   out->Write();
   cout << " WRITE FILE  " << endl;
   out->Close();
   }
//void mass_L0_CCu(TString infile = "out_4.0GeV_CCu_pi3_p4_nMC2000_rebin.root")
 //void mass_L0_CCu(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CCu_pi3_p4_nMC2000_periodI_rebin.root")
   void mass_L0_CCu(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CCu_pi3_p4_nMC2000_periodII_rebin.root")
{

   // ##Read Tree - start
   TFile *file = new TFile(infile, "READ");
   TTree *tree = (TTree *)file->Get("data_inf");

   tree->SetBranchAddress("data_mass", &ms);
   tree->SetBranchAddress("data_pt", &pts);
   tree->SetBranchAddress("data_y", &ys);

   TObjArray *Hlist = new TObjArray();

   for (Int_t i = 0; i < bins_rebin; i++)
   {
      for (Int_t j = 0; j < bins_rebin; j++)
      {
         hname.Form("hdata_rebin%d%d", i + 1, j + 1);
         htitle.Form("%3.2f<y<%3.2f, %3.2f<pt<%3.2f", y_intervals_rebin[i][0], y_intervals_rebin[i][1], pt_intervals_rebin[j][0], pt_intervals_rebin[j][1]);
         sig_data_rebin[i][j] = new TH1F(hname, "", Nbins, min_mass, max_mass);
         sig_data_rebin[i][j]->SetTitle(htitle);
         sig_data_rebin[i][j]->SetMarkerStyle(20);
          sig_data_rebin[i][j]->Sumw2();
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
      // if (pts < 0.1 && pts > 1.05 || ys < 1.2 && ys > 2.1)
      //  continue;
      Int_t i_sig = -999;
      Int_t j_sig = -999;
      if (ys >= 1.2 && ys <= 1.33)
      {
         y_i = 0;
      }
      else if (ys >= 1.33 && ys <= 1.45)
      {
         y_i = 1;
      }
      else if (ys >= 1.45 && ys <= 1.55)
      {
         y_i = 2;
      }
      else if (ys >= 1.55 && ys <= 1.65)
      {
         y_i = 3;
      }
      else if (ys >= 1.65 && ys <= 1.75)
      {
         y_i = 4;
      }
      else if (ys >= 1.75 && ys <= 1.85)
      {
         y_i = 5;
      }
      else if (ys >= 1.85 && ys <= 1.98)
      {
         y_i = 6;
      }
      else if (ys >= 1.98 && ys <= 2.1)
      {
         y_i = 7;
      }
      /*  else
        {
           continue;
        }*/

      if (pts >= 0.1 && pts <= 0.2)
      {
         pt_i = 0;
      }
      else if (pts >= 0.2 && pts <= 0.3)
      {
         pt_i = 1;
      }
      else if (pts >= 0.3 && pts <= 0.4)
      {
         pt_i = 2;
      }
      else if (pts >= 0.4 && pts <= 0.5)
      {
         pt_i = 3;
      }
      else if (pts >= 0.5 && pts <= 0.62)
      {
         pt_i = 4;
      }
      else if (pts >= 0.62 && pts <= 0.75)
      {
         pt_i = 5;
      }
      else if (pts >= 0.75 && pts <= 0.90)
      {
         pt_i = 6;
      }
      else if (pts >= 0.90 && pts <= 1.05)
      {
         pt_i = 7;
      }
      /* else
       {
          continue;
       }*/
      w_cu = 1. / coeff_ccu[pt_i][y_i];
      
       if (w_cu > 100)
         continue;

      //  fill mass spectras in pt&y intervals
      for (Int_t pti = 0; pti < 8; pti++)
      {
         if (pts >= pt_intervals_rebin[pti][0] && pts <= pt_intervals_rebin[pti][1])
         {
            i_sig = pti + 1;
            break;
         }
      }
      for (Int_t yj = 0; yj < 8; yj++)
      {
         if (ys >= y_intervals_rebin[yj][0] && ys <= y_intervals_rebin[yj][1])
         {
            j_sig = yj + 1;

            break;
         }
      }
      if (i_sig > 0 && j_sig > 0)
      {
         sig_data_rebin[j_sig - 1][i_sig - 1]->Fill(ms, w_cu);
         }
 //  }
   
}

Float_t Val[8][8];
Float_t Val_Sum;

for (Int_t k = 0; k < 8; k++)
{
         for (Int_t kk = 0; kk < 8; kk++)
         {
            width_bins_pt[k] = pt_intervals_rebin[k][1] - pt_intervals_rebin[k][0];
            width_bins_y[kk] = y_intervals_rebin[kk][1] - y_intervals_rebin[kk][0];

            Val[k][kk] = width_bins_pt[k] * width_bins_y[kk];
            Val_Sum += Val[k][kk];
           // cout << "  " << width_bins_pt[k] << " ,  " << width_bins_y[kk] << " , " << Val[k][kk] <<endl;
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
            // Hlist->Add(sig_data_rebin_clon[i][j]);

            for (Int_t k = 0; k < sig_data_rebin[i][j]->GetNbinsX(); k++)
            {
               ent = sig_data_rebin[i][j]->GetEntries();

               if(ent < 2)continue;
               // cout << "ent =  " << ent << endl;
               s = sig_data_rebin[i][j]->GetBinContent(k + 1);
               err = sig_data_rebin[i][j]->GetBinError(k + 1);
               if (s == 0 && err == 0)
                  continue;
           //    cout << "  " << s << " , " << err << endl;

               ds = s * Val[i][j] / Val_Sum;
               derr = err * Val[i][j] / Val_Sum;

               if (ds == 0 && derr == 0)
                  continue;
               //cout << "  " << ds << " , " << derr << endl;
               // if(i>=j){
               sig_data_rebin_clon[i][j]->SetBinContent(k + 1, ds);
               sig_data_rebin_clon[i][j]->SetBinError(k + 1, derr);
               //  }
               //}
            }
         }
}
 //###################Procedure extrapolation##################################
for (Int_t ik = 0; ik < bins_rebin; ik++)
{
         for (Int_t jk = 0; jk < bins_rebin; jk++)
         {
            sig_data_rebin_clon_extr[ik][jk] = (TH1F *)sig_data_rebin_clon[ik][jk]->Clone("sig_data_rebin_clon_extr");
         }
}
 ///Double_t wextrp_con[4]={1.493700503 ,1.976310226 ,1.426400599 ,1.024563627};//C+Cu
Double_t wextrp_con[3]={1.4829789,1.193953916,1.042141};//C+Cu
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

             // cout << " getcontent for extrapolate " << s_mass_data[ki][kj] << endl;
             //  cout << "Error for extrapolate " << err_mass_data[ki][kj] << endl;
               
               if (ki == 0)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[0] * s_mass_data[ki][kj];
                  //bin_err = wextrp_con_CCu[0] * err_mass_data[ki][kj];
               }
               else if (ki == 1)
               {
                 biny = wextrp_con[0] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[1] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[1] * err_mass_data[ki][kj];
               }
               else if (ki == 2)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[2] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[2] * err_mass_data[ki][kj];
               } else if (ki == 3)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[2] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[2] * err_mass_data[ki][kj];
               }
                else if (ki == 6)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                  //biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
               } else if (ki == 7)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                  //biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
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
/*for (Int_t ki = 0; ki < bins_rebin; ki++)
{
      for (Int_t kj = 0; kj < bins_rebin; kj++)
      {
      for (Int_t h = 0; h < sig_data_rebin_clon[ki][kj]->GetNbinsX(); h++)
      {
               s_mass_data[ki][kj] = sig_data_rebin[ki][kj]->GetBinContent(h + 1);
               err_mass_data[ki][kj] = sig_data_rebin[ki][kj]->GetBinError(h + 1);

               if(s_mass_data[ki][kj]==0 && err_mass_data[ki][kj]==0) continue;

             // cout << " getcontent for extrapolate " << s_mass_data[ki][kj] << endl;
             //  cout << "Error for extrapolate " << err_mass_data[ki][kj] << endl;
               
               if (ki == 0)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[0] * s_mass_data[ki][kj];
                  //bin_err = wextrp_con_CCu[0] * err_mass_data[ki][kj];
               }
               else if (ki == 1)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[1] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[1] * err_mass_data[ki][kj];
               }
               else if (ki == 2)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[2] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[2] * err_mass_data[ki][kj];
               }
                else if (ki == 7)
               {
                 biny = wextrp_con[3] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[3] * err_mass_data[ki][kj];
                  //biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
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
}*/
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
  // hist_sum_by_pt->Sumw2();
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
  // hist_sum_by_y->Sumw2();
   Hlist->Add(hist_sum_by_y);
 //=================Sum Hist суммируем строчки//
    TH1F *hpt_b1 = (TH1F *)summ_4pt_b1->Clone("hpt_b1");
   TH1F *hpt_b2 = (TH1F *)summ_4pt_b2->Clone("hpt_b2");
   TH1F *hpt_b3 = (TH1F *)summ_4pt_b3->Clone("hpt_b3");
   TH1F *hpt_b4 = (TH1F *)summ_4pt_b4->Clone("hpt_b4");

   Double_t s_pt,s_pt2,s_pt3,s_pt4;
   Double_t err_pt,err_pt2,err_pt3,err_pt4;
   Double_t bin1,bin2,bin3,bin4;
   Double_t bin_err1,bin_err2,bin_err3,bin_err4;

   for (Int_t h = 0; h < summ_4pt_b1->GetNbinsX(); h++)
   {
      s_pt = summ_4pt_b1->GetBinContent(h + 1);
      err_pt = summ_4pt_b1->GetBinError(h + 1);

      bin1 =0.58*s_pt;
      bin_err1=0.58*err_pt;

      hpt_b1->SetBinContent(h + 1, bin1);
      hpt_b1->SetBinError(h + 1, bin_err1);
   }
     for (Int_t h1 = 0; h1 < summ_4pt_b2->GetNbinsX(); h1++)
   {
      s_pt2 = summ_4pt_b2->GetBinContent(h1 + 1);
      err_pt2 = summ_4pt_b2->GetBinError(h1 + 1);

      bin2 =0.59*s_pt2;
      bin_err2=0.59*err_pt2;

      hpt_b2->SetBinContent(h1 + 1, bin2);
      hpt_b2->SetBinError(h1 + 1, bin_err2);
   }

    for (Int_t h2 = 0; h2< summ_4pt_b3->GetNbinsX(); h2++)
   {
      s_pt3 = summ_4pt_b3->GetBinContent(h2 + 1);
      err_pt3 = summ_4pt_b3->GetBinError(h2 + 1);

      bin3 =0.58*s_pt3;
      bin_err3=0.58*err_pt3;

      hpt_b3->SetBinContent(h2 + 1, bin3);
      hpt_b3->SetBinError(h2 + 1, bin_err3);
   }
    for (Int_t h3 = 0; h3 < summ_4pt_b4->GetNbinsX(); h3++)
   {
      s_pt4 = summ_4pt_b4->GetBinContent(h3 + 1);
      err_pt4 = summ_4pt_b4->GetBinError(h3 + 1);

      bin4 =0.56*s_pt4;
      bin_err4=0.56*err_pt4;

      hpt_b4->SetBinContent(h3 + 1, bin4);
      hpt_b4->SetBinError(h3 + 1, bin_err4);
   }

   Hlist->Add(hpt_b1);
   Hlist->Add(hpt_b2);
   Hlist->Add(hpt_b3);
   Hlist->Add(hpt_b4);

   //TFile *out = new TFile("outFile_full_8bins_part_CCu_all.root", "RECREATE");
    //TFile *out = new TFile("outFile_full_8bins_part_CCu_full_refit.root", "RECREATE");
    //TFile *out = new TFile("outFile_full_8bins_part_CCu_full_refit_perI.root", "RECREATE");
    TFile *out = new TFile("outFile_full_8bins_part_CCu_full_refit_perII.root", "RECREATE");
   //TFile *out = new TFile("outFile_full_8bins_part_CCu_periodI.root", "RECREATE");
   //TFile *out = new TFile("outFile_full_8bins_part_CCu_periodII.root", "RECREATE");
   Hlist->Write();
   out->Write();
   cout << " WRITE FILE  " << endl;
   out->Close();
   }

  void mass_L0_CAl(TString infile = "out_4.0GeV_CAl_pi3_p4_nMC2000_rebin.root")
  //  void mass_L0_CAl(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CAl_pi3_p4_nMC2000_periodI_rebin.root")
 // void mass_L0_CAl(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CAl_pi3_p4_nMC2000_periodII_rebin.root")
{

   // ##Read Tree - start
   TFile *file = new TFile(infile, "READ");
   TTree *tree = (TTree *)file->Get("data_inf");
   tree->SetBranchAddress("data_mass", &ms);
   tree->SetBranchAddress("data_pt", &pts);
   tree->SetBranchAddress("data_y", &ys);

   Int_t pt_i, y_i;
   // коэфефициенты
  // Float_t w_cc, w_cu, w_cal, w_cpb;
   TObjArray *Hlist = new TObjArray();
   // ###################### rebin
   for (Int_t i = 0; i < bins_rebin; i++)
   {
      for (Int_t j = 0; j < bins_rebin; j++)
      {
         hname.Form("hdata_rebin%d%d", i + 1, j + 1);
         htitle.Form("%3.2f<y<%3.2f, %3.2f<pt<%3.2f", y_intervals_rebin[i][0], y_intervals_rebin[i][1], pt_intervals_rebin[j][0], pt_intervals_rebin[j][1]);
         sig_data_rebin[i][j] = new TH1F(hname, "", Nbins, min_mass, max_mass);
         sig_data_rebin[i][j]->SetTitle(htitle);
         sig_data_rebin[i][j]->SetMarkerStyle(20);
         sig_data_rebin[i][j]->Sumw2();
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
      // if (pts < 0.1 && pts > 1.05 || ys < 1.2 && ys > 2.1)
      //  continue;

      Int_t i_sig = -999;
      Int_t j_sig = -999;
      if (ys >= 1.2 && ys <= 1.33)
      {
         y_i = 0;
      }
      else if (ys >= 1.33 && ys <= 1.45)
      {
         y_i = 1;
      }
      else if (ys >= 1.45 && ys <= 1.55)
      {
         y_i = 2;
      }
      else if (ys >= 1.55 && ys <= 1.65)
      {
         y_i = 3;
      }
      else if (ys >= 1.65 && ys <= 1.75)
      {
         y_i = 4;
      }
      else if (ys >= 1.75 && ys <= 1.85)
      {
         y_i = 5;
      }
      else if (ys >= 1.85 && ys <= 1.98)
      {
         y_i = 6;
      }
      else if (ys >= 1.98 && ys <= 2.1)
      {
         y_i = 7;
      }
      /*  else
        {
           continue;
        }*/

      if (pts >= 0.1 && pts <= 0.2)
      {
         pt_i = 0;
      }
      else if (pts >= 0.2 && pts <= 0.3)
      {
         pt_i = 1;
      }
      else if (pts >= 0.3 && pts <= 0.4)
      {
         pt_i = 2;
      }
      else if (pts >= 0.4 && pts <= 0.5)
      {
         pt_i = 3;
      }
      else if (pts >= 0.5 && pts <= 0.62)
      {
         pt_i = 4;
      }
      else if (pts >= 0.62 && pts <= 0.75)
      {
         pt_i = 5;
      }
      else if (pts >= 0.75 && pts <= 0.90)
      {
         pt_i = 6;
      }
      else if (pts >= 0.90 && pts <= 1.05)
      {
         pt_i = 7;
      }
      /* else
       {
          continue;
       }*/
       w_cal = 1./coeff_cal[pt_i][y_i];
       if (w_cal > 100)continue;
      //  fill mass spectras in pt&y intervals
      for (Int_t pti = 0; pti < 8; pti++)
      {
         if (pts >= pt_intervals_rebin[pti][0] && pts <= pt_intervals_rebin[pti][1])
         {
            i_sig = pti + 1;
            break;
         }
      }
      for (Int_t yj = 0; yj < 8; yj++)
      {
         if (ys >= y_intervals_rebin[yj][0] && ys <= y_intervals_rebin[yj][1])
         {
            j_sig = yj + 1;

            break;
         }
      }
      if (i_sig > 0 && j_sig > 0)
      {
         sig_data_rebin[j_sig - 1][i_sig - 1]->Fill(ms, w_cal);
      }
 //  }
   
}

Float_t Val[8][8];
Float_t Val_Sum;

for (Int_t k = 0; k < 8; k++)
{
         for (Int_t kk = 0; kk < 8; kk++)
         {
            width_bins_pt[k] = pt_intervals_rebin[k][1] - pt_intervals_rebin[k][0];
            width_bins_y[kk] = y_intervals_rebin[kk][1] - y_intervals_rebin[kk][0];

            Val[k][kk] = width_bins_pt[k] * width_bins_y[kk];
            Val_Sum += Val[k][kk];
           // cout << "  " << width_bins_pt[k] << " ,  " << width_bins_y[kk] << " , " << Val[k][kk] <<endl;
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
            // Hlist->Add(sig_data_rebin_clon[i][j]);

            for (Int_t k = 0; k < sig_data_rebin[i][j]->GetNbinsX(); k++)
            {
               ent = sig_data_rebin[i][j]->GetEntries();

               if(ent < 2)continue;
               // cout << "ent =  " << ent << endl;
               s = sig_data_rebin[i][j]->GetBinContent(k + 1);
               err = sig_data_rebin[i][j]->GetBinError(k + 1);
               if (s == 0 && err == 0)
                  continue;
           //    cout << "  " << s << " , " << err << endl;

               ds = s * Val[i][j] / Val_Sum;
               derr = err * Val[i][j] / Val_Sum;

               if (ds == 0 && derr == 0)
                  continue;
               //cout << "  " << ds << " , " << derr << endl;
               // if(i>=j){
               sig_data_rebin_clon[i][j]->SetBinContent(k + 1, ds);
               sig_data_rebin_clon[i][j]->SetBinError(k + 1, derr);
               //  }
               //}
            }
         }
}
 //###################Procedure extrapolation##################################
//wextrp_con[4];

for (Int_t ik = 0; ik < bins_rebin; ik++)
{
         for (Int_t jk = 0; jk < bins_rebin; jk++)
         {

            sig_data_rebin_clon_extr[ik][jk] = (TH1F *)sig_data_rebin_clon[ik][jk]->Clone("sig_data_rebin_clon_extr");
            // sig_data_rebin_clon_extr[ik][jk]->Reset();
            //    sig_data_rebin_clon[i][j]= (TH1F*) sig_data_rebin[i][j]->Clone("sig_data_rebin_clon");
            // Hlist->Add(sig_data_rebin_clon_extr[ik][jk]);
         }
}
//Double_t wextrp_con[5] = {2.037187652, 1.48024508, 1.061358536, 1.049894024, 1.020058267}; // C+Al
Double_t wextrp_con[3] = {1.726205429, 1.055985585, 1.007500823}; // C+Al
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

             // cout << " getcontent for extrapolate " << s_mass_data[ki][kj] << endl;
             //  cout << "Error for extrapolate " << err_mass_data[ki][kj] << endl;
               
               if (ki == 0)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[0] * s_mass_data[ki][kj];
                  //bin_err = wextrp_con_CCu[0] * err_mass_data[ki][kj];
               }
               else if (ki == 1)
               {
                 biny = wextrp_con[0] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[1] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[1] * err_mass_data[ki][kj];
               }
               else if (ki == 2)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[2] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[2] * err_mass_data[ki][kj];
               }
               else if (ki == 3)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[3] * s_mass_data[ki][kj];
                  //bin_err = wextrp_con_CCu[3] * err_mass_data[ki][kj];
               }
               else if (ki == 6)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
               }
               else if (ki == 7)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
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
   // hist_sum_by_pt->Sumw2();
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
  // hist_sum_by_y->Sumw2();
   Hlist->Add(hist_sum_by_y);

  // TFile *out = new TFile("outFile_full_8bins_part_CAl_periodII.root", "RECREATE");
   //TFile *out = new TFile("outFile_full_8bins_part_CAl_periodI.root", "RECREATE");
   TFile *out = new TFile("outFile_full_8bins_part_CAl_all.root", "RECREATE");
   Hlist->Write();
   out->Write();
   cout << " WRITE FILE  " << endl;
   out->Close();
   }

   void mass_L0_CPb(TString infile = "/home/alish/Data_C+T_r6/KA/out_4.0GeV_CPb_pi3_p4_nMC2000_periodI_rebin.root")
   //void mass_L0_CPb(TString infile = "out_4.0GeV_CPb_pi3_p4_nMC2000_rebin.root")
{

   // ##Read Tree - start
   TFile *file = new TFile(infile, "READ");
   TTree *tree = (TTree *)file->Get("data_inf");

   Double_t pts, ys, ms;

   tree->SetBranchAddress("data_mass", &ms);
   tree->SetBranchAddress("data_pt", &pts);
   tree->SetBranchAddress("data_y", &ys);

   TObjArray *Hlist = new TObjArray();
   // ###################### rebin
   for (Int_t i = 0; i < bins_rebin; i++)
   {
      for (Int_t j = 0; j < bins_rebin; j++)
      {
         hname.Form("hdata_rebin%d%d", i + 1, j + 1);
         htitle.Form("%3.2f<y<%3.2f, %3.2f<pt<%3.2f", y_intervals_rebin[i][0], y_intervals_rebin[i][1], pt_intervals_rebin[j][0], pt_intervals_rebin[j][1]);
         sig_data_rebin[i][j] = new TH1F(hname, "", Nbins, min_mass, max_mass);
         sig_data_rebin[i][j]->SetTitle(htitle);
         sig_data_rebin[i][j]->SetMarkerStyle(20);
         sig_data_rebin[i][j]->Sumw2();
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
      // if (pts < 0.1 && pts > 1.05 || ys < 1.2 && ys > 2.1)
      //  continue;

      Int_t i_sig = -999;
      Int_t j_sig = -999;
      if (ys >= 1.2 && ys <= 1.33)
      {
         y_i = 0;
      }
      else if (ys >= 1.33 && ys <= 1.45)
      {
         y_i = 1;
      }
      else if (ys >= 1.45 && ys <= 1.55)
      {
         y_i = 2;
      }
      else if (ys >= 1.55 && ys <= 1.65)
      {
         y_i = 3;
      }
      else if (ys >= 1.65 && ys <= 1.75)
      {
         y_i = 4;
      }
      else if (ys >= 1.75 && ys <= 1.85)
      {
         y_i = 5;
      }
      else if (ys >= 1.85 && ys <= 1.98)
      {
         y_i = 6;
      }
      else if (ys >= 1.98 && ys <= 2.1)
      {
         y_i = 7;
      }
      /*  else
        {
           continue;
        }*/

      if (pts >= 0.1 && pts <= 0.2)
      {
         pt_i = 0;
      }
      else if (pts >= 0.2 && pts <= 0.3)
      {
         pt_i = 1;
      }
      else if (pts >= 0.3 && pts <= 0.4)
      {
         pt_i = 2;
      }
      else if (pts >= 0.4 && pts <= 0.5)
      {
         pt_i = 3;
      }
      else if (pts >= 0.5 && pts <= 0.62)
      {
         pt_i = 4;
      }
      else if (pts >= 0.62 && pts <= 0.75)
      {
         pt_i = 5;
      }
      else if (pts >= 0.75 && pts <= 0.90)
      {
         pt_i = 6;
      }
      else if (pts >= 0.90 && pts <= 1.05)
      {
         pt_i = 7;
      }
      /* else
       {
          continue;
       }*/
       w_cpb = 1./ coeff_cpb[pt_i][y_i];

       if (w_cpb > 100)continue;
      //  fill mass spectras in pt&y intervals
      for (Int_t pti = 0; pti < 8; pti++)
      {
         if (pts >= pt_intervals_rebin[pti][0] && pts <= pt_intervals_rebin[pti][1])
         {
            i_sig = pti + 1;
            break;
         }
      }
      for (Int_t yj = 0; yj < 8; yj++)
      {
         if (ys >= y_intervals_rebin[yj][0] && ys <= y_intervals_rebin[yj][1])
         {
            j_sig = yj + 1;

            break;
         }
      }
      if (i_sig > 0 && j_sig > 0)
      {
         sig_data_rebin[j_sig - 1][i_sig - 1]->Fill(ms, w_cpb);
      }
 //  }
   
}
Float_t Val[8][8];
Float_t Val_Sum;

for (Int_t k = 0; k < 8; k++)
{
         for (Int_t kk = 0; kk < 8; kk++)
         {
           width_bins_pt[k] = pt_intervals_rebin[k][1] - pt_intervals_rebin[k][0];
           width_bins_y[kk] = y_intervals_rebin[kk][1] - y_intervals_rebin[kk][0];

           Val[k][kk] = width_bins_pt[k] * width_bins_y[kk];
           Val_Sum += Val[k][kk];
           // cout << "  " << width_bins_pt[k] << " ,  " << width_bins_y[kk] << " , " << Val[k][kk] <<endl;
         }
}
// cout << " sum = " <<  Val_Sum << endl;
//  TO

//cout << " sum = " <<  Val_Sum << endl;
// TObjArray *Hlist = new TObjArray();
Double_t s, err, ds, derr;
Double_t ent;
for (Int_t i = 0; i < bins_rebin; i++)
{
         for (Int_t j = 0; j < bins_rebin; j++)
         {
            sig_data_rebin_clon[i][j] = (TH1F *)sig_data_rebin[i][j]->Clone("sig_data_rebin_clon");
            // Hlist->Add(sig_data_rebin_clon[i][j]);

            for (Int_t k = 0; k < sig_data_rebin[i][j]->GetNbinsX(); k++)
            {
               ent = sig_data_rebin[i][j]->GetEntries();

               if(ent < 2)continue;
               // cout << "ent =  " << ent << endl;
               s = sig_data_rebin[i][j]->GetBinContent(k + 1);
               err = sig_data_rebin[i][j]->GetBinError(k + 1);
               if (s == 0 && err == 0)
                  continue;
           //    cout << "  " << s << " , " << err << endl;

               ds = s * Val[i][j] / Val_Sum;
               derr = err * Val[i][j] / Val_Sum;

               if (ds == 0 && derr == 0)
                  continue;
               //cout << "  " << ds << " , " << derr << endl;
               // if(i>=j){
               sig_data_rebin_clon[i][j]->SetBinContent(k + 1, ds);
               sig_data_rebin_clon[i][j]->SetBinError(k + 1, derr);
               //  }
               //}
            }
         }
}
 //###################Procedure extrapolation##################################
for (Int_t ik = 0; ik < bins_rebin; ik++)
{
         for (Int_t jk = 0; jk < bins_rebin; jk++)
         {
            sig_data_rebin_clon_extr[ik][jk] = (TH1F *)sig_data_rebin_clon[ik][jk]->Clone("sig_data_rebin_clon_extr");
         }
} 
Double_t wextrp_con[3] = {2.390230448,1.130666295,1.028777987}; // C+Pb
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

             // cout << " getcontent for extrapolate " << s_mass_data[ki][kj] << endl;
             //  cout << "Error for extrapolate " << err_mass_data[ki][kj] << endl;
               
               if (ki == 0)
               {
                  biny = wextrp_con[0] * s_mass_data[ki][kj];
                  bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[0] * s_mass_data[ki][kj];
                  //bin_err = wextrp_con_CCu[0] * err_mass_data[ki][kj];
               }
               else if (ki == 1)
               {
                 biny = wextrp_con[0] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[0] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[1] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[1] * err_mass_data[ki][kj];
               }
               else if (ki == 2)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[2] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[2] * err_mass_data[ki][kj];
               }
               else if (ki == 3)
               {
                 biny = wextrp_con[1] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[1] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[3] * s_mass_data[ki][kj];
                  //bin_err = wextrp_con_CCu[3] * err_mass_data[ki][kj];
               }
               else if (ki == 6)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
               }
               else if (ki == 7)
               {
                 biny = wextrp_con[2] * s_mass_data[ki][kj];
                 bin_err = wextrp_con[2] * err_mass_data[ki][kj];
                 // biny = wextrp_con_CCu[4] * s_mass_data[ki][kj];
                 // bin_err = wextrp_con_CCu[4] * err_mass_data[ki][kj];
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
   // hist_sum_by_pt->Sumw2();
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
   //hist_sum_by_y->Sumw2();
   Hlist->Add(hist_sum_by_y);

   TFile *out = new TFile("outFile_full_8bins_part_CPb_periodI.root", "RECREATE");
   Hlist->Write();
   out->Write();
   cout << " WRITE FILE  " << endl;
   out->Close();
}