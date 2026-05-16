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


void tsachi1()
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

// 参数顺序：
// [0] = Nor
// [1] = q
// [2] = T

double NorMin = 1.0;      // 自己改范围
double NorMax = 50.0;     // 自己改范围
int nScan = 1000;

double best_q = 0;
double best_Nor = 0;
double best_T = 0;
double best_chi2 = 1e99;

int N = gr->GetN();

// 初值
double q_init = 1.11;
double T_init = 0.10;

TGraph *gChi2Nor = new TGraph();
int ip = 0;

for (int i = 0; i <= nScan; i++) {

    double Norval = NorMin + (NorMax - NorMin) * i / nScan;

    // -------------------------
    // 设置参数初值
    // -------------------------
    f1->SetParameter(0, Norval);
    f1->SetParameter(1, q_init);
    f1->SetParameter(2, T_init);

    // -------------------------
    // 设置物理范围
    // -------------------------
    f1->SetParLimits(0, 0.0, 1e6);      // Nor > 0
    f1->SetParLimits(1, 1.01, 1.25);    // q 范围
    f1->SetParLimits(2, 0.001, 1.0);    // T > 0

    // -------------------------
    // 固定 Nor，拟合 q 和 T
    // -------------------------
    f1->FixParameter(0, Norval);
    f1->ReleaseParameter(1);
    f1->ReleaseParameter(2);

    TFitResultPtr r = gr->Fit(f1, "S Q R");

    // 拟合失败就跳过
    if (int(r) != 0) continue;

    double q = f1->GetParameter(1);
    double T = f1->GetParameter(2);

    // 非物理结果跳过
    if (q <= 1.0 || T <= 0) continue;

    // -------------------------
    // 用统计误差重新计算 chi2
    // -------------------------

	
	double chi2 = f1->GetChisquare();
	
	
	gChi2Nor->SetPoint(ip, Norval, chi2);
    	ip++;
    
    
    // -------------------------
    // 找最小统计 chi2
    // -------------------------
    if (chi2 < best_chi2) {

        best_chi2 = chi2;
        best_Nor = Norval;
        best_q = q;
        best_T = T;
    }

    // 用本次结果作为下一次初值，提高稳定性
    q_init = q;
    T_init = T;
}

// 自由度：N - 自由参数数目
// 这里 Nor 固定，q 和 T 自由，所以自由参数数目是 2
int ndf = N - 2;

cout << "Best Nor = " << best_Nor << endl;
cout << "Best q = " << best_q << endl;
cout << "Best T = " << best_T << endl;
cout << "Best chi2 = " << best_chi2 << endl;
cout << "Best chi2/ndf = " << best_chi2 / ndf << endl;
		


	
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


TCanvas *cChi2 = new TCanvas("cChi2", "chi2 vs Nor", 550, 500);

gChi2Nor->SetTitle("#chi^{2}_{stat} vs Nor;Nor;#chi^{2}_{stat}");
gChi2Nor->SetMarkerStyle(20);
gChi2Nor->SetMarkerSize(0.8);
gChi2Nor->SetLineWidth(2);

gChi2Nor->Draw("ALP");

// 画出最小值位置
TLine *lineBest = new TLine(best_Nor, gChi2Nor->GetYaxis()->GetXmin(),
                            best_Nor, gChi2Nor->GetYaxis()->GetXmax());
lineBest->SetLineStyle(2);
lineBest->SetLineWidth(2);
lineBest->Draw("same");


cChi2->SaveAs("chi2_vs_Nor.png");
		
// chi2	
	
	Int_t N1 = gr->GetN();
	
	Int_t ndf1 = N1 - parN;

	double chi2 = 0;
	double chi2st = 0;
	double chi2sy = 0;
	
	for(int i = 1; i < gr->GetN() + 1; i++)
	{
		
		chi2 = chi2 + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / (pow(e1->GetBinError(i), 2) + pow(e2->GetBinError(i), 2));
	
		chi2sy = chi2sy + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / pow(e2->GetBinError(i), 2);
	
		chi2st = chi2st + pow(f1->Eval(gr->GetPointX(i-1)) / 1 - gr->GetPointY(i-1) / 1, 2) / pow(e1->GetBinError(i), 2);
	
		
	
	}
	
	
	
	cout << "N = " << N1 << endl;
	cout << "chi2 = " << chi2 << endl;
	cout << "chi2sy = " << chi2sy << endl;
	cout << "chi2st = " << chi2st << endl;
	cout << "ndf = " << ndf1 <<endl;
}


