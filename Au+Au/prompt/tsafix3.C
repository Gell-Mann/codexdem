#include <TMath.h>
#include <math.h>

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


void tsafix3()
{
	TFile *f = TFile::Open("HEPData-ins2061074-v1-Figure_15.root");
	
	TDirectoryFile *td = (TDirectoryFile*)f->Get("Figure 15");
	
	TGraphAsymmErrors *gr = (TGraphAsymmErrors*)td->Get("Graph1D_y3");
	
	TH1F *e1 = (TH1F*)td->Get("Hist1D_y3_e1");
	TH1F *e2 = (TH1F*)td->Get("Hist1D_y3_e2");
	
	
	TCanvas *c1 = new TCanvas("fit1","fit1",960,0,550,500);
	
	double xMin = 0.;
	double xMax = 10.;
	int parN = 3;   
        
        
        TF1 *f1 =new TF1("f1", func1, xMin, xMax, parN);

	f1->SetParNames("Nor", "q", "T");

	f1->FixParameter(0, 32.202);
     	f1->FixParameter(1, 1.130);
	f1->FixParameter(2, 0.059);
	
	// 65.891
	// 0.003	
	// 0.017
	
	c1->SetTickx();
	c1->SetTicky();
	c1->SetLogy();
		
	gr->SetLineColor(kCyan-2);
	gr->SetLineWidth(2);
	gr->SetMarkerStyle(20);
	gr->SetMarkerColor(4);
	gr->SetMarkerSize(1);
	
	gr->GetXaxis()->SetLimits(0, 12);
	
	gr->SetMinimum(pow(10, -8));
	gr->SetMaximum(pow(10, 2));

		
	gStyle->SetTitleY(0.96);
	gr->SetTitle("");
	
	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitle("p_{T} (GeV)");
	gr->GetXaxis()->SetTitleOffset(1.3);
	
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetLabelSize(0.03);
	gr->GetYaxis()->SetTitle("1/N 1/(2#pip_{T}) d^{2}N/(dydp_{T}) (GeV/c)^{-2}");
	gr->GetYaxis()->SetTitleOffset(1.3);
	
	gr->Draw("AP");
	gr->Fit("f1");
		
		
// chi2	
	
	Int_t N = gr->GetN();
	
	Int_t ndf = N - parN;

	double chi2 = 0;
	
	for(int i = 1; i < gr->GetN() + 1; i++)
	{
		
		chi2 = chi2 + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / (pow(e1->GetBinError(i), 2) + pow(e2->GetBinError(i), 2));
	
	}
	
	
	
	cout << "N = " << N << endl;
	cout << "chi2 = " << chi2 << endl;
	cout << "ndf = " << ndf <<endl;
}
