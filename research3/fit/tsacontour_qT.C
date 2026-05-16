#include <TMath.h>
#include <math.h>

using namespace TMath;

Double_t func1(Double_t *x, Double_t *par)
{
    Double_t a = 0;

    Double_t b = TMath::Sqrt(TMath::Power(x[0], 2) + TMath::Power(a, 2));

    Double_t z1 = TMath::Power(1 + (par[1] - 1) * b / par[2],
                                par[1] / (1 - par[1]));

    Double_t z2 = 2 * par[0] * b / TMath::Power(2 * TMath::Pi(), 3);

    return z1 * z2 * TMath::Power(10, 4);
}


void tsacontour_qT()
{
    TFile *f = TFile::Open("HEPData-ins2061074-v1-Figure_15.root");

    TDirectoryFile *td = (TDirectoryFile*)f->Get("Figure 15");

    TGraphAsymmErrors *gr = (TGraphAsymmErrors*)td->Get("Graph1D_y1");

    TH1F *e1 = (TH1F*)td->Get("Hist1D_y1_e1");
    TH1F *e2 = (TH1F*)td->Get("Hist1D_y1_e2");

    double xMin = 0.;
    double xMax = 10.;
    int parN = 3;

    TF1 *f1 = new TF1("f1_qT_contour", func1, xMin, xMax, parN);

    // 参数顺序：
    // [0] = Nor
    // [1] = q
    // [2] = T

    int N = gr->GetN();

    // -----------------------------
    // q 和 T 的扫描范围
    // -----------------------------
    double qMin = 1.08;
    double qMax = 1.15;
    int nq = 120;

    double TMin = 0.05;
    double TMax = 0.20;
    int nT = 120;

    // -----------------------------
    // 最优值保存
    // -----------------------------
    double best_q = 0;
    double best_T = 0;
    double best_Nor = 0;
    double best_chi2st = 1e99;

    double Nor_init = 10.0;

    // 二维 chi2 图
    TH2D *hChi2 = new TH2D("hChi2",
                           "#chi^{2}_{stat}(q,T);q;T",
                           nq, qMin, qMax,
                           nT, TMin, TMax);

    // -----------------------------
    // 二维扫描 q 和 T
    // 每个点固定 q 和 T，只拟合 Nor
    // -----------------------------
    for (int iq = 0; iq < nq; iq++) {

        double qval = qMin + (qMax - qMin) * iq / (nq - 1);

        for (int iT = 0; iT < nT; iT++) {

            double Tval = TMin + (TMax - TMin) * iT / (nT - 1);

            f1->SetParameter(0, Nor_init);
            f1->SetParameter(1, qval);
            f1->SetParameter(2, Tval);

            f1->SetParLimits(0, 0.0, 1e6);
            f1->SetParLimits(1, 1.01, 1.25);
            f1->SetParLimits(2, 0.001, 1.0);

            // 固定 q 和 T，只拟合 Nor
            f1->ReleaseParameter(0);
            f1->FixParameter(1, qval);
            f1->FixParameter(2, Tval);

            TFitResultPtr r = gr->Fit(f1, "S Q R N");

            if (int(r) != 0) {
                hChi2->SetBinContent(iq + 1, iT + 1, 1e9);
                continue;
            }

            double Nor = f1->GetParameter(0);

            if (Nor <= 0) {
                hChi2->SetBinContent(iq + 1, iT + 1, 1e9);
                continue;
            }

            // -----------------------------
            // 用统计误差重新计算 chi2
            // -----------------------------
            double chi2st = 0.0;

            for (int j = 0; j < N; j++) {

                double x, y;
                gr->GetPoint(j, x, y);

                double yfit = f1->Eval(x);

                double ey_stat = e1->GetBinError(j + 1);

                if (ey_stat <= 0) continue;

                chi2st += pow((y - yfit) / ey_stat, 2);
            }

            hChi2->SetBinContent(iq + 1, iT + 1, chi2st);

            if (chi2st < best_chi2st) {
                best_chi2st = chi2st;
                best_q = qval;
                best_T = Tval;
                best_Nor = Nor;
            }

            Nor_init = Nor;
        }
    }

    int ndf = N - 1;   // q 和 T 固定，只拟合 Nor，一个自由参数

    cout << "Best q = " << best_q << endl;
    cout << "Best T = " << best_T << endl;
    cout << "Best Nor = " << best_Nor << endl;
    cout << "Best chi2_stat = " << best_chi2st << endl;
    cout << "Best chi2_stat/ndf = " << best_chi2st / ndf << endl;
    cout << "ndf = " << ndf << endl;


    // -----------------------------
    // 画二维彩色图
    // -----------------------------
    TCanvas *c2 = new TCanvas("c2", "chi2 q T contour", 550, 500);

    gStyle->SetOptStat(0);

    TH2D *hDraw = (TH2D*)hChi2->Clone("hDraw");
hDraw->SetTitle("#Delta#chi^{2}_{stat}(q,T);q;T");

for (int ix = 1; ix <= hDraw->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= hDraw->GetNbinsY(); iy++) {

        double v = hDraw->GetBinContent(ix, iy);

        if (v <= 0 || v > best_chi2st + 500) {
            hDraw->SetBinContent(ix, iy, 500);
        } else {
            hDraw->SetBinContent(ix, iy, v - best_chi2st);
        }
    }
}

hDraw->SetMinimum(0);
hDraw->SetMaximum(20);
    
    
hDraw->Draw("COLZ");

TMarker *mBest = new TMarker(best_q, best_T, 29);
mBest->SetMarkerColor(kRed);
mBest->SetMarkerSize(2.0);
mBest->Draw("same");

c2->Update();
c2->SaveAs("chi2_qT_colz.png");


    // -----------------------------
    // 画等高线图
    // -----------------------------
    TCanvas *c3 = new TCanvas("c3", "chi2 q T contour lines", 550, 500);

    double levels[6];

    levels[0] = best_chi2st + 1.0;
    levels[1] = best_chi2st + 2.30;
    levels[2] = best_chi2st + 4.0;
    levels[3] = best_chi2st + 6.18;
    levels[4] = best_chi2st + 9.0;
    levels[5] = best_chi2st + 11.83;

    hChi2->SetContour(6, levels);

    hChi2->Draw("CONT1");

    TMarker *mBest2 = new TMarker(best_q, best_T, 29);
    mBest2->SetMarkerColor(kRed);
    mBest2->SetMarkerSize(2.0);
    mBest2->Draw("same");

    c3->SaveAs("chi2_qT_contour.png");


    // -----------------------------
    // 画最佳拟合谱
    // -----------------------------
    TCanvas *c1 = new TCanvas("fit_best", "best fit", 960, 0, 550, 500);

    c1->SetTickx();
    c1->SetTicky();
    c1->SetLogy();

    gr->SetLineColor(kCyan - 2);
    gr->SetLineWidth(2);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(1);

    gr->GetXaxis()->SetLimits(0, 12);

    gr->SetMinimum(pow(10, -8));
    gr->SetMaximum(pow(10, 2));

    gr->SetTitle("");

    gr->GetXaxis()->CenterTitle();
    gr->GetXaxis()->SetTitle("p_{T} (GeV)");
    gr->GetXaxis()->SetTitleOffset(1.3);

    gr->GetYaxis()->CenterTitle();
    gr->GetYaxis()->SetLabelSize(0.03);
    gr->GetYaxis()->SetTitle("1/N 1/(2#pip_{T}) d^{2}N/(dydp_{T}) (GeV/c)^{-2}");
    gr->GetYaxis()->SetTitleOffset(1.3);

    gr->Draw("AP");

    f1->SetParameter(0, best_Nor);
    f1->SetParameter(1, best_q);
    f1->SetParameter(2, best_T);

    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->Draw("same");

    c1->SaveAs("best_fit_qT.png");
}
