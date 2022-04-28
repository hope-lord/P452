# include <iostream>
# include <cmath>
# include <cstring>

# include <TApplication.h>
# include <TCanvas.h>
# include <TMultiGraph.h>
# include <TGraph.h>
# include <TLegend.h>

using std::string;

void plot(string f1, string f2, string f3){
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TMultiGraph *mg = new TMultiGraph();
	TGraph *g1 = new TGraph(f1.c_str(), "%lg %lg");
	TGraph *g2 = new TGraph(f2.c_str(), "%lg %lg");
	TGraph *g3 = new TGraph(f3.c_str(), "%lg %lg");
	g1->SetMarkerColor(kRed);
	g2->SetMarkerColor(kBlue);
	g3->SetMarkerColor(kGreen);
	g1->SetMarkerStyle(29);
	g2->SetMarkerStyle(5);
	g3->SetMarkerStyle(22);
	g1->SetLineColor(kRed);
	g2->SetLineColor(kBlue);
	g3->SetLineColor(kGreen);
	g1->SetLineWidth(1);
	g2->SetLineWidth(1);
	g3->SetLineWidth(1);
	mg->Add(g1);
	mg->Add(g2);
	mg->Add(g3);
	mg->SetTitle("Wave Function;x;#psi(x)");
	TLegend leg(0.1, 0.7, 0.3, 0.9);
	leg.AddEntry(g1, "Guess 1", "lp");
	leg.AddEntry(g2, "Solution", "lp");
	leg.AddEntry(g3, "Guess 2", "lp");
	mg->Draw("ACP");
	leg.DrawClone("Same");
}

int main(int argc,char*argv[]){
	if (argc !=4){
		printf( "Usage: ./plotq2 <file1> <file2> <file3>\n" );
		return 0;
	}
	string f1 = argv[1];
	string f2 = argv[2];
	string f3 = argv[3];
	TApplication *app = new TApplication("app)", &argc, argv);
	plot(f1,f2,f3);
	app->Run();
	return 0;
}


