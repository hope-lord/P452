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
  TGraph *gr1 = new TGraph("q2_Jacobi_iter.txt","%lg %lg");
  TGraph *gr2 = new TGraph("q2_Gauss_Seidal_iter.txt","%lg %lg");
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
  c1->SaveAs("q2_comparision.png");

  TCanvas *c2 = new TCanvas("c2","Residual");
  TMultiGraph *mg2 = new TMultiGraph();
  c2->SetLogy(1);
	TGraph *g = new TGraph("q3_residual.txt","%lg %lg");
	mg2->SetTitle("For Q3 Conjugate Gradient Method;Iteration;Residual in log scale");
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  mg2->Add(g);
  mg2->Draw("ACP");
  c2->SaveAs("q3_res.png");
  gSystem->ProcessEvents();
  return 0;
}
