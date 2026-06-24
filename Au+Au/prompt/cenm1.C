#include "TMath.h"
#include "TFile.h"
#include <Math/Functor.h>
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "Fit/Fitter.h"
#include "iostream"

using namespace TMath;

Double_t func1(Double_t *x, Double_t *par)
{
	Double_t a = 0;
	
	Double_t z1;
	Double_t b;
	Double_t z2;

	b = TMath::Sqrt(TMath::Power(x[0], 2) + TMath::Power(a, 2));

	z1 = TMath::Power(1 + (par[1] - 1) * b / par[2], par[1] / (1 - par[1]));

	z2 = 2 * par[0] * b / TMath::Power(2 * TMath::Pi(), 3);

	return z1 * z2 * TMath::Power(10, 4);

	
}




void cenm1()
{
	TFile *f = TFile::Open("HEPData-ins2061074-v1-Figure_15.root");
	
	TDirectoryFile *td = (TDirectoryFile*)f->Get("Figure 15");
	
	TGraphAsymmErrors *gr1 = (TGraphAsymmErrors*)td->Get("Graph1D_y1");
	
	TGraphAsymmErrors *gr2 = (TGraphAsymmErrors*)td->Get("Graph1D_y2");
	
	TGraphAsymmErrors *gr3 = (TGraphAsymmErrors*)td->Get("Graph1D_y3");
	
	TGraphAsymmErrors *gr4 = (TGraphAsymmErrors*)td->Get("Graph1D_y4");	
	
	
	
	TCanvas *c1 = new TCanvas("fit1","fit1",960,0,550,500);
	

	double xMin = 0.;
	double xMax = 12.;
	int parN = 3;   
        TF1 *f1 =new TF1("f1", func1, xMin, xMax, parN);
        
        f1->SetLineColor(2);
        
        f1->SetParNames("Nor", "q", "T");
        
        
	f1->FixParameter(0, 30*8.363);
     	f1->SetParameter(1, 1.110);
	f1->FixParameter(2, 0.107);
     	        
	for(int i = 0; i < gr1->GetN(); i++)
	{
		double x, y;
		gr1->GetPoint(i, x, y);
		
		gr1->SetPoint(i, x, 30*y);
		gr1->SetPointEYlow(i,30*gr1->GetErrorYlow(i));
		gr1->SetPointEYhigh(i, 30*gr1->GetErrorYhigh(i));
	}        

	
	c1->SetTickx();
	c1->SetTicky();
	c1->SetLogy();
	
	gr1->GetXaxis()->SetLimits(0, 11);
	gr1->SetMinimum(pow(10, -9));
	gr1->SetMaximum(pow(10, 3));
	
	gr1->SetLineColor(kCyan-2);
	gr1->SetLineWidth(2);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(4);
	gr1->SetMarkerSize(1);
	
	
	
	gStyle->SetTitleY(0.96);
	gStyle->SetLineScalePS(2);
	
	gStyle->SetLineStyleString(10, "15 15");


	
	gr1->SetTitle("");
	
	
	gr1->GetXaxis()->CenterTitle();
	gr1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr1->GetXaxis()->SetTitleOffset(1.3);
	
	gr1->GetYaxis()->CenterTitle();
	gr1->GetYaxis()->SetLabelSize(0.03);
	gr1->GetYaxis()->SetTitle("1/N 1/(2#pip_{T}) d^{2}N/(dydp_{T}) ((GeV/c)^{-2})");
	gr1->GetYaxis()->SetTitleOffset(1.3);
	
	gr1->Draw("AP");
	gr1->Fit("f1");
	f1->Draw("same");

		
	
// second curve
	
	TF1 *f2 =new TF1("f2", func1, xMin, xMax, parN);
        
        f2->SetLineColor(2);
       
        
	f2->FixParameter(0, 10*21.143);
     	f2->SetParameter(1, 1.120);
	f2->FixParameter(2, 0.078);

	for(int i = 0; i < gr2->GetN(); i++)
	{
		double x, y;
		gr2->GetPoint(i, x, y);
		
		gr2->SetPoint(i, x, 10*y);
		gr2->SetPointEYlow(i,10*gr2->GetErrorYlow(i));
		gr2->SetPointEYhigh(i, 10*gr2->GetErrorYhigh(i));
	}  
     	     	   
	
	gr2->SetLineColor(kCyan-2);
	gr2->SetLineWidth(2);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(1);
	gr2->SetMarkerSize(1);
	
	gr2->Draw("P same");
	gr2->Fit("f2");
	f2->Draw("same");
	
		
	
// third curve
	
	TF1 *f3 =new TF1("f3", func1, xMin, xMax, parN);
        
        f3->SetLineColor(2);
        
        
	f3->FixParameter(0, 5*32.202);
     	f3->FixParameter(1, 1.130);
	f3->FixParameter(2, 0.059);

	for(int i = 0; i < gr3->GetN(); i++)
	{
		double x, y;
		gr3->GetPoint(i, x, y);
		
		gr3->SetPoint(i, x, 5*y);
		gr3->SetPointEYlow(i,5*gr3->GetErrorYlow(i));
		gr3->SetPointEYhigh(i, 5*gr3->GetErrorYhigh(i));
	}       	               
	
	gr3->SetLineColor(kCyan-2);
	gr3->SetLineWidth(2);
	gr3->SetMarkerStyle(22);
	gr3->SetMarkerColor(4);
	gr3->SetMarkerSize(1);
	
	gr3->Draw("P same");
	gr3->Fit("f3");
	f3->Draw("same");
			
	
// fourth curve
	
	TF1 *f4 =new TF1("f4", func1, xMin, xMax, parN);
        
        f4->SetLineColor(2);
        
	f4->FixParameter(0, 40.185);
     	f4->SetParameter(1, 1.142);
	f4->FixParameter(2, 0.040);        

	for(int i = 0; i < gr4->GetN(); i++)
	{
		double x, y;
		gr4->GetPoint(i, x, y);
		
		gr4->SetPoint(i, x, 1*y);
		gr4->SetPointEYlow(i,1*gr4->GetErrorYlow(i));
		gr4->SetPointEYhigh(i, 1*gr4->GetErrorYhigh(i));
	}  
            	     
	
	gr4->SetLineColor(kCyan-2);
	gr4->SetLineWidth(2);
	gr4->SetMarkerStyle(23);
	gr4->SetMarkerColor(1);
	gr4->SetMarkerSize(1);
	
	gr4->Draw("P same");
	gr4->Fit("f4");
	f4->Draw("same");
	
			
	
	
// fit function style
   	
   	TLine *line1 = new TLine(0.20, 0.78, 0.25, 0.78); line1->SetNDC();
   	
   	TLatex *tex1 = new TLatex; tex1->SetNDC();
   	tex1->SetTextSize(0.03);
   	tex1->SetTextFont(42);
   	tex1->DrawLatex(0.27,  0.77, "Tsallis");
   	
   	
   	line1->SetLineStyle(1);
   	line1->SetLineColor(2);
   	line1->SetLineWidth(2);
   	line1->Draw("same");
   	

// data entry  	
   	
   	TLatex *tex20 = new TLatex; tex20->SetNDC();
   	tex20->SetTextSize(0.03);
   	tex20->SetTextFont(132);
   	tex20->DrawLatex(0.2,  0.84, "direct-photon 200 GeV Au+Au, PHENIX");
   	
   	
   	TLine *line21 = new TLine(0.65, 0.85, 0.70, 0.85); line21->SetNDC();
   	line21->SetLineColor(kCyan-2);
	line21->SetLineWidth(2);
   	line21->Draw();
   	
   	TMarker *ma = new TMarker((0.65 + 0.70) / 2, 0.85, 20); ma->SetNDC();
   	
   	ma->SetMarkerColor(4);
	ma->SetMarkerSize(1);
	ma->Draw("same");
	
	TLatex *tex21 = new TLatex; tex21->SetNDC();
   	tex21->SetTextSize(0.03);
   	tex21->SetTextFont(42);
   	tex21->DrawLatex(0.71,  0.84, "0#minus20\%, #times30");
   	
   	TLine *line22 = new TLine(0.65, 0.80, 0.70, 0.80); line22->SetNDC();
   	line22->SetLineColor(kCyan-2);
	line22->SetLineWidth(2);
   	line22->Draw();
   	
   	TMarker *ma2 = new TMarker((0.65 + 0.70) / 2, 0.80, 21); ma2->SetNDC();
   	
   	ma2->SetMarkerColor(1);
	ma2->SetMarkerSize(1);
	ma2->Draw("same");
	
	TLatex *tex22 = new TLatex; tex22->SetNDC();
   	tex22->SetTextSize(0.03);
   	tex22->SetTextFont(42);
   	tex22->DrawLatex(0.71,  0.79, "20#minus40\%, #times10");
   	
   	TLine *line23 = new TLine(0.65, 0.75, 0.70, 0.75); line23->SetNDC();
   	line23->SetLineColor(kCyan-2);
	line23->SetLineWidth(2);
   	line23->Draw();
   	
   	TMarker *ma3 = new TMarker((0.65 + 0.70) / 2, 0.75, 22); ma3->SetNDC();
   	
   	ma3->SetMarkerColor(4);
	ma3->SetMarkerSize(1);
	ma3->Draw("same");
	
	TLatex *tex23 = new TLatex; tex23->SetNDC();
   	tex23->SetTextSize(0.03);
   	tex23->SetTextFont(42);
   	tex23->DrawLatex(0.71,  0.74, "40#minus60\%, #times5");
   	
   	TLine *line24 = new TLine(0.65, 0.70, 0.70, 0.70); line24->SetNDC();
   	line24->SetLineColor(kCyan-2);
	line24->SetLineWidth(2);
   	line24->Draw();
   	
   	TMarker *ma4= new TMarker((0.65 + 0.70) / 2, 0.70, 23); ma4->SetNDC();
   	
   	ma4->SetMarkerColor(1);
	ma4->SetMarkerSize(1);
	ma4->Draw("same");
	
	TLatex *tex24 = new TLatex; tex24->SetNDC();
   	tex24->SetTextSize(0.03);
   	tex24->SetTextFont(42);
   	tex24->DrawLatex(0.71,  0.69, "60#minus93\%, #times1");
   	
/*   	TLine *line25 = new TLine(0.60, 0.65, 0.65, 0.65); line25->SetNDC();
   	line25->SetLineColor(kCyan-2);
	line25->SetLineWidth(2);
   	line25->Draw();
   	
   	TMarker *ma5= new TMarker((0.60 + 0.65) / 2, 0.65, 29); ma5->SetNDC();
   	
   	ma5->SetMarkerColor(4);
	ma5->SetMarkerSize(1.2);
	ma5->Draw("same");
	
	TLatex *tex25 = new TLatex; tex25->SetNDC();
   	tex25->SetTextSize(0.03);
   	tex25->SetTextFont(42);
   	tex25->DrawLatex(0.66,  0.64, "60#minus92\%, #times1");
*/
   	   	
/*   	TLatex *tex26 = new TLatex; tex26->SetNDC();
   	tex26->SetTextSize(0.03);
   	tex26->SetTextFont(42);
   	tex26->DrawLatex(0.66,  0.59, "|y| < 0.5");
*/
  	
	TLatex *tex27 = new TLatex; tex27->SetNDC();
   	tex27->SetTextSize(0.04);
   	tex27->SetTextFont(42);
   	tex27->DrawLatex(0.15,  0.15, "(a)");
   	    		
	
}

