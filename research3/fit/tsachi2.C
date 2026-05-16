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


void tsachi2()
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

double qMin = 1.08;
double qMax = 1.14;
int nScan = 1000;

double best_q = 0;
double best_Nor = 0;
double best_T = 0;
double best_chi2 = 1e99;

int N = gr->GetN();

// 初值
double Nor_init = 10.0;
double T_init   = 0.10;

TGraph *gChi2q = new TGraph();
int ip = 0;

for (int i = 0; i <= nScan; i++) {

    double qval = qMin + (qMax - qMin) * i / nScan;

    // -------------------------
    // 设置参数初值
    // -------------------------
    f1->SetParameter(0, Nor_init);
    f1->SetParameter(1, qval);
    f1->SetParameter(2, T_init);

    // -------------------------
    // 设置物理范围
    // -------------------------
    f1->SetParLimits(0, 0.0, 1e6);      // Nor > 0
    f1->SetParLimits(1, 1.01, 1.25);    // q 范围
    f1->SetParLimits(2, 0.001, 1.0);    // T > 0

    // -------------------------
    // 固定 q，拟合 Nor 和 T
    // -------------------------
    f1->ReleaseParameter(0);
    f1->FixParameter(1, qval);
    f1->ReleaseParameter(2);

    TFitResultPtr r = gr->Fit(f1, "S Q R");
    gr->Draw("AP");

    if (int(r) != 0) continue;

    double Nor = f1->GetParameter(0);
    double T   = f1->GetParameter(2);

    if (Nor <= 0 || T <= 0) continue;

    // -------------------------
    // 用统计误差重新计算 chi2
    // -------------------------
	
	double chi2 = f1->GetChisquare();

    gChi2q->SetPoint(ip, qval, chi2);
    ip++;

    // -------------------------
    // 找最小统计 chi2
    // -------------------------
    if (chi2 < best_chi2) {

        best_chi2 = chi2;
        best_q = qval;
        best_Nor = Nor;
        best_T = T;
    }

    // 用本次结果作为下一次初值
    Nor_init = Nor;
    T_init   = T;
}

// 自由度：q 固定，Nor 和 T 自由
int ndf = N - 2;

cout << "Best q = " << best_q << endl;
cout << "Best Nor = " << best_Nor << endl;
cout << "Best T = " << best_T << endl;
cout << "Best chi2 = " << best_chi2 << endl;
cout << "Best chi2/ndf = " << best_chi2 / ndf << endl;
		
		f1->SetParameter(0, best_Nor);
		f1->SetParameter(1, best_q);
		f1->SetParameter(2, best_T);

	
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
	
	//gr->Draw("AP");


TCanvas *cChi2 = new TCanvas("cChi2", "chi2 vs q", 800, 600);

gChi2q->SetTitle("#chi^{2}_{stat} vs q;q;#chi^{2}_{stat}");
gChi2q->SetMarkerStyle(20);
gChi2q->SetMarkerSize(0.8);
gChi2q->SetLineWidth(2);

gChi2q->Draw("ALP");

// 画出最小值位置
TLine *lineBest = new TLine(best_q,
                            gChi2q->GetYaxis()->GetXmin(),
                            best_q,
                            gChi2q->GetYaxis()->GetXmax());
lineBest->SetLineStyle(2);
lineBest->SetLineWidth(2);
lineBest->Draw("same");

cChi2->SaveAs("chi2_vs_q.png");
		
// chi2	
	
	Int_t N1 = gr->GetN();
	
	Int_t ndf1 = N1 - parN;

	double chi2 = 0;
	double chi2st = 0;
	double chi2sy = 0;
	
	for(int i = 0; i < gr->GetN(); i++)
	{
		

		
		//chi2 = chi2 + pow(f1->Eval(gr->GetPointX(i-1)) - gr->GetPointY(i-1), 2) / (pow(e1->GetBinError(i), 2) + pow(e2->GetBinError(i), 2));
	
		//chi2sy = chi2sy + pow(f1->Eval(gr->GetPointX(i-1)) - gr->GetPointY(i-1), 2) / pow(e2->GetBinError(i), 2);
	
		chi2st = chi2st + pow(f1->Eval(gr->GetPointX(i)) - gr->GetPointY(i), 2) / pow(e1->GetBinError(i+1), 2);
	
		
	
	}
	
	
	
	cout << "N = " << N1 << endl;
	cout << "chi2 = " << chi2 << endl;
	cout << "chi2sy = " << chi2sy << endl;
	cout << "chi2st = " << chi2st << endl;
	cout << "ndf = " << ndf1 <<endl;
}


