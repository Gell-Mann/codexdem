#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "iostream"

void sca2()
{

	
	double x1[] = {10, 30, 50, 76.5};
	double y1[] = {0.107, 0.078, 0.059, 0.040};
        
        double ex1[] = {0., 0., 0., 0., 0.};
        double ey1[] = {0.017, 0.017, 0.017, 0.019};
        
        int numpoints = 4;
        
        
        TGraphErrors *graph1 = new TGraphErrors(numpoints, x1, y1, ex1, ey1);
	
	
	TCanvas *c1 = new TCanvas("sca","sca",960,570,500,450);
	
	gStyle->SetTitleY(0.96);
	graph1->SetTitle("");
	
	graph1->GetXaxis()->CenterTitle();
	graph1->GetXaxis()->SetTitle("Centrality");
	graph1->GetXaxis()->SetTitleOffset(1.3);
	
	graph1->GetYaxis()->CenterTitle();
	graph1->GetYaxis()->SetLabelSize(0.03);
	graph1->GetYaxis()->SetTitle("T_{eff}");
	graph1->GetYaxis()->SetTitleOffset(1.48);
	
	//graph1->SetMinimum(0.3);
	//graph1->SetMaximum(0.4);
	
	graph1->SetMarkerColor(4);
	graph1->SetMarkerStyle(20);
	graph1->SetMarkerSize(1.3);
	graph1->Draw("AP");
        
        c1->SetTickx();
        c1->SetTicky();
        


// data entry           
        
   	TLatex *tex1 = new TLatex; tex1->SetNDC();
   	tex1->SetTextSize(0.04);
   	tex1->SetTextFont(42);
   	tex1->DrawLatex(0.17,  0.84, "200 GeV Au+Au");
   	
     	TLatex *tex2 = new TLatex; tex2->SetNDC();
   	tex2->SetTextSize(0.04);
   	tex2->SetTextFont(42);
   	tex2->DrawLatex(0.15,  0.15, "(c)");
        
}
