#include <TMath.h>
#include <TError.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <math.h>
#include <iostream>
#include <memory>

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


void tsa4()
{
	constexpr const char *inputPath = "HEPData-ins2061074-v1-Figure_15.root";
	constexpr const char *directoryName = "Figure 15";
	constexpr const char *graphName = "Graph1D_y4";
	constexpr const char *statErrorName = "Hist1D_y4_e1";
	constexpr const char *sysErrorName = "Hist1D_y4_e2";
	
	const int oldErrorIgnoreLevel = gErrorIgnoreLevel;
	gErrorIgnoreLevel = kFatal;
	
	std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
	if (!input || input->IsZombie()) {
		gErrorIgnoreLevel = oldErrorIgnoreLevel;
		std::cerr << "Error: cannot open " << inputPath << '\n';
		return;
	}
	
	auto *td = dynamic_cast<TDirectoryFile*>(input->Get(directoryName));
	if (!td) {
		gErrorIgnoreLevel = oldErrorIgnoreLevel;
		std::cerr << "Error: cannot find directory \"" << directoryName << "\"\n";
		return;
	}
	
	auto *inputGraph = dynamic_cast<TGraphAsymmErrors*>(td->Get(graphName));
	auto *inputE1 = dynamic_cast<TH1F*>(td->Get(statErrorName));
	auto *inputE2 = dynamic_cast<TH1F*>(td->Get(sysErrorName));
	if (!inputGraph || !inputE1 || !inputE2) {
		gErrorIgnoreLevel = oldErrorIgnoreLevel;
		std::cerr << "Error: cannot find graph or error histograms\n";
		return;
	}
	
	auto *gr = dynamic_cast<TGraphAsymmErrors*>(inputGraph->Clone("directPhotonSpectrum_y4"));
	auto *e1 = dynamic_cast<TH1F*>(inputE1->Clone("directPhotonSpectrum_y4_stat"));
	auto *e2 = dynamic_cast<TH1F*>(inputE2->Clone("directPhotonSpectrum_y4_sys"));
	if (!gr || !e1 || !e2) {
		gErrorIgnoreLevel = oldErrorIgnoreLevel;
		std::cerr << "Error: cannot clone graph or error histograms\n";
		return;
	}
	e1->SetDirectory(nullptr);
	e2->SetDirectory(nullptr);

	input.reset();
	gErrorIgnoreLevel = oldErrorIgnoreLevel;
	
	
	TCanvas *c1 = new TCanvas("fit1","fit1",960,0,550,500);
	
	double xMin = 0.;
	double xMax = 10.;
	int parN = 3;   
        
        
        TF1 *f1 =new TF1("f1", func1, xMin, xMax, parN);

	f1->SetParNames("Nor", "q", "T");


	f1->SetParameter(0, 0.312712);
	f1->SetParLimits(0, 0., 2.5);
     	f1->SetParameter(1, 1.142);
	f1->SetParLimits(1, 1.001, 1.5);
	f1->FixParameter(2, 0.073);
		
	
	c1->SetTickx();
	c1->SetTicky();
	c1->SetLogy();
		
	gr->SetLineColor(kCyan-2);
	gr->SetLineWidth(2);
	gr->SetMarkerStyle(20);
	gr->SetMarkerColor(4);
	gr->SetMarkerSize(1);
	
	gr->GetXaxis()->SetLimits(0, 12);
	
	gr->SetMinimum(pow(10, -9));
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
	gErrorIgnoreLevel = kFatal;
	gr->Fit("f1", "R");
		
		
// chi2	
	
	Int_t N = 0;
	double chi2 = 0;
	double chi2st = 0;
	
	for(int i = 1; i < gr->GetN() + 1; i++)
	{
		double x = gr->GetPointX(i-1);
		if (x < xMin || x > xMax) continue;
		N++;
		chi2 = chi2 + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / (pow(e1->GetBinError(i), 2) + pow(e2->GetBinError(i), 2));
		chi2st = chi2st + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / pow(e1->GetBinError(i), 2);
	
	}
	
	Int_t ndf = N - parN;
	
	
	
	cout << "N = " << N << endl;
	cout << "Nor = " << f1->GetParameter(0) << " +/- " << f1->GetParError(0) << endl;
	cout << "q = " << f1->GetParameter(1) << " +/- " << f1->GetParError(1) << endl;
	cout << "T = " << f1->GetParameter(2) << " +/- " << f1->GetParError(2) << endl;
	cout << "chi2 = " << chi2 << endl;
	cout << "chi2st = " << chi2st << endl;
	cout << "ndf = " << ndf <<endl;
}
