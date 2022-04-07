#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TGraph.h>
using namespace std;

int main(int argc,char ** argv){
TApplication app("app",&argc,argv);
  TCanvas *c1 = new TCanvas("c1","Comparision");
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Coparision between Gauss-Seidel and Jacobi Mehthod;Iteration No.;Residual in log scale");
  c1->SetLogy(1);
  TGraph *gr1 = new TGraph("Jacobi_iter.txt","%lg %lg");
  TGraph *gr2 = new TGraph("Gauss_Seidal_iter.txt","%lg %lg");
  gr1->SetMarkerStyle(20);
  gr2->SetMarkerStyle(21);
  gr1->SetMarkerColor(kRed);
  gr2->SetMarkerColor(kBlue);
  gr1->SetLineColor(kRed);
  gr2->SetLineColor(kBlue);
  mg -> Add(gr1);
  mg -> Add(gr2);
  /* gr1->GetXaxis()->SetTitle("Iteration"); */
  /* gr1->GetYaxis()->SetTitle("Error"); */
  TLegend leg (0.9,0.8,0.6,0.9);
  leg.AddEntry(gr1,"Jacobi Method","lp");
  leg.AddEntry(gr2,"Gauss-Seidal Method","lp");
  mg->Draw("ACP");
  leg.DrawClone("same");
  c1->SaveAs("comparision.pdf");
  gSystem->ProcessEvents();
  return 0;
}
