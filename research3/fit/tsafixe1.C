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


void tsafixe1()
{
	TFile *f = TFile::Open("HEPData-ins2061074-v1-Figure_15.root");
	
	TDirectoryFile *td = (TDirectoryFile*)f->Get("Figure 15");
	
	TGraphAsymmErrors *gr = (TGraphAsymmErrors*)td->Get("Graph1D_y1");
	
	TH1F *e1 = (TH1F*)td->Get("Hist1D_y1_e1");
	TH1F *e2 = (TH1F*)td->Get("Hist1D_y1_e2");
	
	
	TCanvas *c1 = new TCanvas("fit1","fit1",960,0,550,500);
	
	double xMin = 0.;
	double xMax = 10.;
	int parN = 3;   
        
        
        TF1 *f1 =new TF1("f1", func1, xMin, xMax, parN);

// 假设参数顺序：
// [0] = Nor
// [1] = q
// [2] = T

// -------------------------
// 先定义初值
// -------------------------
double Nor0 = 8.363;
double q0   = 1.1109;
double T0   = 0.1066;

// =====================================================
// 1. 固定 Nor，拟合 q 和 T，求 q、T 的 MINOS 误差
// =====================================================
f1->SetParameter(0, Nor0);
f1->SetParameter(1, q0);
f1->SetParameter(2, T0);

f1->FixParameter(0, Nor0);      // 固定 Nor
f1->ReleaseParameter(1);        // 放开 q
f1->ReleaseParameter(2);        // 放开 T

TFitResultPtr r_qT = gr->Fit(f1, "S E");

double q_best  = r_qT->Parameter(1);
double q_err_l = fabs(r_qT->LowerError(1));
double q_err_h = r_qT->UpperError(1);

double T_best  = r_qT->Parameter(2);
double T_err_l = fabs(r_qT->LowerError(2));
double T_err_h = r_qT->UpperError(2);


// =====================================================
// 2. 固定 q 和 T，拟合 Nor，求 Nor 的 MINOS 误差
// =====================================================
f1->SetParameter(0, Nor0);
f1->SetParameter(1, q_best);
f1->SetParameter(2, T_best);

f1->ReleaseParameter(0);        // 放开 Nor
f1->FixParameter(1, q_best);    // 固定 q
f1->FixParameter(2, T_best);    // 固定 T

TFitResultPtr r_Nor = gr->Fit(f1, "S E");

double Nor_best  = r_Nor->Parameter(0);
double Nor_err_l = fabs(r_Nor->LowerError(0));
double Nor_err_h = r_Nor->UpperError(0);


// =====================================================
// 3. 输出结果
// =====================================================
cout << "Nor = " << Nor_best
     << " +" << Nor_err_h
     << " -" << Nor_err_l << endl;

cout << "q = " << q_best
     << " +" << q_err_h
     << " -" << q_err_l << endl;

cout << "T = " << T_best
     << " +" << T_err_h
     << " -" << T_err_l << endl;
		
	
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
/*	TFitResultPtr r = gr->Fit("f1", "S E");
	
	double q      = r->Parameter(1);
	double q_err_low  = r->LowerError(1);
	double q_err_high = r->UpperError(1);

	double T      = r->Parameter(2);
	double T_err_low  = r->LowerError(2);
	double T_err_high = r->UpperError(2);

	cout << "q = " << q 
     	<< " +" << q_err_high 
     	<< " " << q_err_low << endl;

	cout << "T = " << T 
     	<< " +" << T_err_high 
     	<< " " << T_err_low << endl;	
*/		
// chi2	
	
	Int_t N = gr->GetN();
	
	Int_t ndf = N - parN;

	double chi2 = 0;
	double chi2st = 0;
	double chi2sy = 0;
	
	for(int i = 1; i < gr->GetN() + 1; i++)
	{
		
		chi2 = chi2 + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / (pow(e1->GetBinError(i), 2) + pow(e2->GetBinError(i), 2));
	
		chi2sy = chi2sy + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / pow(e2->GetBinError(i), 2);
	
		chi2st = chi2st + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / pow(e1->GetBinError(i), 2);
	
		
	
	}
	
	
	
	cout << "N = " << N << endl;
	cout << "chi2 = " << chi2 << endl;
	cout << "chi2sy = " << chi2sy << endl;
	cout << "chi2st = " << chi2st << endl;
	cout << "ndf = " << ndf <<endl;
}



