# include <iostream>

# include <TApplication.h>
# include <TCanvas.h>
# include <TGraph2D.h>

void plot(){
  TCanvas *c = new TCanvas("c", "3D Plot", 800, 600);
  TGraph2D *g = new TGraph2D("q3_out.txt");
  g->SetTitle("#frac{#partial^{2}#phi}{#partialx^{2}} + #frac{#partial^{2}#phi}{#partialy^{2}} = 0, #phi(x,0)=1, #phi(x,1)=#phi(0,y)=#phi(1,y)=0;x;y;#phi(x,y)");
  g->Draw("surf1");
}

int main(int argc,char*argv[]){
	TApplication app("app",&argc,argv);
	plot();
	app.Run();
	return 0;
}
