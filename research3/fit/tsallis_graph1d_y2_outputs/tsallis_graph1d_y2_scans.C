#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TError.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

Double_t TsallisFuncAll(Double_t *x, Double_t *par)
{
    const Double_t a = 0.0;
    const Double_t b = TMath::Sqrt(TMath::Power(x[0], 2) + TMath::Power(a, 2));
    const Double_t z1 = TMath::Power(1 + (par[1] - 1) * b / par[2],
                                     par[1] / (1 - par[1]));
    const Double_t z2 = 2 * par[0] * b / TMath::Power(2 * TMath::Pi(), 3);

    return z1 * z2 * TMath::Power(10, 4);
}

struct ScanPoint {
    double fixedValue;
    double Nor;
    double q;
    double T;
    double chi2;
};

struct ScanResult {
    string name;
    string parameterName;
    int fixedIndex;
    vector<ScanPoint> points;
    int bestIndex;
    double bestNor;
    double bestQ;
    double bestT;
    double bestChi2;
    double errLow;
    double errHigh;
    double parErrLow[3];
    double parErrHigh[3];
    bool minosOk;
    int ndf;
};

struct ContourResult {
    double bestNor;
    double bestQ;
    double bestT;
    double bestChi2;
    int ndf;
};

struct PairProfileResult {
    string name;
    string xName;
    string yName;
    int xIndex;
    int yIndex;
    int freeIndex;
    double bestPars[3];
    double parErrLow[3];
    double parErrHigh[3];
    double errorDeltaChi2;
    double bestChi2;
    int ndf;
};

struct PairProfilePoint {
    double pars[3];
    double chi2;
    bool valid;

    PairProfilePoint() : chi2(0.0), valid(false)
    {
        for (int ipar = 0; ipar < 3; ++ipar) {
            pars[ipar] = 0.0;
        }
    }
};

const double kNorMin = 0.30;
const double kNorMax = 40.0;
const double kQMin = 1.08;
const double kQMax = 1.14;
const double kTMin = 0.05;
const double kTMax = 0.20;
const int kScan1D = 1000;
const int kNQContour = 120;
const int kNTContour = 120;
const double kBadChi2 = 1.0e30;
const double kPairErrorDeltaChi2 = 1.0;
const string kOutputDir = ".";

string OutputPath(const string &fileName)
{
    return kOutputDir + "/" + fileName;
}

void ConfigureFunction(TF1 *f)
{
    f->SetParName(0, "Nor");
    f->SetParName(1, "q");
    f->SetParName(2, "T");
    f->SetParLimits(0, 0.0, 1.0e6);
    f->SetParLimits(1, 1.01, 1.25);
    f->SetParLimits(2, 0.001, 1.0);
}

int CountStatPoints(TGraphAsymmErrors *gr, TH1F *eStat)
{
    int n = 0;
    for (int i = 0; i < gr->GetN(); ++i) {
        if (eStat->GetBinError(i + 1) > 0) {
            ++n;
        }
    }
    return n;
}

TGraphErrors *BuildStatFitGraph(TGraphAsymmErrors *gr, TH1F *eStat)
{
    TGraphErrors *gFit = new TGraphErrors();
    gFit->SetName("gFit_stat_errors");

    int ip = 0;
    for (int i = 0; i < gr->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        gr->GetPoint(i, x, y);

        const double ey = eStat->GetBinError(i + 1);
        if (ey <= 0) {
            continue;
        }

        gFit->SetPoint(ip, x, y);
        gFit->SetPointError(ip, 0.0, ey);
        ++ip;
    }

    return gFit;
}

double Chi2Stat(TGraphAsymmErrors *gr, TH1F *eStat, TF1 *f)
{
    double chi2 = 0.0;

    for (int i = 0; i < gr->GetN(); ++i) {
        double x = 0.0;
        double y = 0.0;
        gr->GetPoint(i, x, y);

        const double ey = eStat->GetBinError(i + 1);
        if (ey <= 0) {
            continue;
        }

        const double diff = y - f->Eval(x);
        chi2 += diff * diff / (ey * ey);
    }

    return chi2;
}

bool FitAtFixedValue(TGraphErrors *gFit,
                     TGraphAsymmErrors *gr,
                     TH1F *eStat,
                     TF1 *f,
                     int fixedIndex,
                     double fixedValue,
                     double &initNor,
                     double &initQ,
                     double &initT,
                     ScanPoint &point)
{
    f->SetParameters(initNor, initQ, initT);
    f->SetParameter(fixedIndex, fixedValue);
    ConfigureFunction(f);

    for (int ipar = 0; ipar < 3; ++ipar) {
        if (ipar == fixedIndex) {
            f->FixParameter(ipar, fixedValue);
        } else {
            f->ReleaseParameter(ipar);
        }
    }

    TFitResultPtr r = gFit->Fit(f, "SQNR");
    if (int(r) != 0) {
        return false;
    }

    const double Nor = f->GetParameter(0);
    const double q = f->GetParameter(1);
    const double T = f->GetParameter(2);

    if (Nor <= 0 || q <= 1.0 || T <= 0) {
        return false;
    }

    const double chi2 = Chi2Stat(gr, eStat, f);

    point.fixedValue = fixedValue;
    point.Nor = Nor;
    point.q = q;
    point.T = T;
    point.chi2 = chi2;

    initNor = Nor;
    initQ = q;
    initT = T;

    return true;
}

void InitializeParameterErrors(ScanResult &result)
{
    for (int i = 0; i < 3; ++i) {
        result.parErrLow[i] = numeric_limits<double>::quiet_NaN();
        result.parErrHigh[i] = numeric_limits<double>::quiet_NaN();
    }
}

void InitializePairProfileErrors(PairProfileResult &result)
{
    result.errorDeltaChi2 = kPairErrorDeltaChi2;

    for (int i = 0; i < 3; ++i) {
        result.parErrLow[i] = numeric_limits<double>::quiet_NaN();
        result.parErrHigh[i] = numeric_limits<double>::quiet_NaN();
    }
}

double InterpolateCrossing(const ScanPoint &a, const ScanPoint &b, double target)
{
    const double dy = b.chi2 - a.chi2;
    if (dy == 0.0) {
        return numeric_limits<double>::quiet_NaN();
    }

    return a.fixedValue + (target - a.chi2) * (b.fixedValue - a.fixedValue) / dy;
}

void UpdatePairParameterRanges(const double pars[3],
                               double minPars[3],
                               double maxPars[3],
                               bool &hasRange)
{
    hasRange = true;

    for (int ipar = 0; ipar < 3; ++ipar) {
        minPars[ipar] = std::min(minPars[ipar], pars[ipar]);
        maxPars[ipar] = std::max(maxPars[ipar], pars[ipar]);
    }
}

void UpdatePairBoundaryRanges(const PairProfilePoint &a,
                              const PairProfilePoint &b,
                              double target,
                              double minPars[3],
                              double maxPars[3],
                              bool &hasRange)
{
    if (!a.valid || !b.valid) {
        return;
    }

    const bool aInside = a.chi2 <= target;
    const bool bInside = b.chi2 <= target;
    if (aInside == bInside) {
        return;
    }

    const double dChi2 = b.chi2 - a.chi2;
    if (dChi2 == 0.0) {
        return;
    }

    const double fraction = (target - a.chi2) / dChi2;
    if (fraction < 0.0 || fraction > 1.0) {
        return;
    }

    double pars[3];
    for (int ipar = 0; ipar < 3; ++ipar) {
        pars[ipar] = a.pars[ipar] + fraction * (b.pars[ipar] - a.pars[ipar]);
    }

    UpdatePairParameterRanges(pars, minPars, maxPars, hasRange);
}

void ComputePairProfileErrors(PairProfileResult &result,
                              const vector<vector<PairProfilePoint> > &grid)
{
    double minPars[3] = {numeric_limits<double>::infinity(),
                         numeric_limits<double>::infinity(),
                         numeric_limits<double>::infinity()};
    double maxPars[3] = {-numeric_limits<double>::infinity(),
                         -numeric_limits<double>::infinity(),
                         -numeric_limits<double>::infinity()};
    bool hasRange = false;

    if (!std::isfinite(result.bestChi2)) {
        return;
    }

    const double target = result.bestChi2 + result.errorDeltaChi2;

    for (int ix = 0; ix < (int)grid.size(); ++ix) {
        for (int iy = 0; iy < (int)grid[ix].size(); ++iy) {
            const PairProfilePoint &point = grid[ix][iy];

            if (point.valid && point.chi2 <= target) {
                UpdatePairParameterRanges(point.pars, minPars, maxPars, hasRange);
            }

            if (ix + 1 < (int)grid.size()) {
                UpdatePairBoundaryRanges(point, grid[ix + 1][iy], target,
                                         minPars, maxPars, hasRange);
            }

            if (iy + 1 < (int)grid[ix].size()) {
                UpdatePairBoundaryRanges(point, grid[ix][iy + 1], target,
                                         minPars, maxPars, hasRange);
            }
        }
    }

    if (!hasRange) {
        return;
    }

    for (int ipar = 0; ipar < 3; ++ipar) {
        result.parErrLow[ipar] = std::max(0.0, result.bestPars[ipar] - minPars[ipar]);
        result.parErrHigh[ipar] = std::max(0.0, maxPars[ipar] - result.bestPars[ipar]);
    }
}

void ComputeOneSigmaErrors(ScanResult &result)
{
    result.errLow = numeric_limits<double>::quiet_NaN();
    result.errHigh = numeric_limits<double>::quiet_NaN();

    if (result.points.empty() || result.bestIndex < 0) {
        return;
    }

    const double target = result.bestChi2 + 1.0;

    for (int i = result.bestIndex; i > 0; --i) {
        const ScanPoint &right = result.points[i];
        const ScanPoint &left = result.points[i - 1];
        const double yRight = right.chi2 - target;
        const double yLeft = left.chi2 - target;

        if (yRight == 0.0) {
            result.errLow = result.points[result.bestIndex].fixedValue - right.fixedValue;
            break;
        }

        if (yLeft == 0.0 || yLeft * yRight < 0.0) {
            const double crossing = InterpolateCrossing(left, right, target);
            result.errLow = result.points[result.bestIndex].fixedValue - crossing;
            break;
        }
    }

    for (int i = result.bestIndex; i + 1 < (int)result.points.size(); ++i) {
        const ScanPoint &left = result.points[i];
        const ScanPoint &right = result.points[i + 1];
        const double yLeft = left.chi2 - target;
        const double yRight = right.chi2 - target;

        if (yLeft == 0.0) {
            result.errHigh = left.fixedValue - result.points[result.bestIndex].fixedValue;
            break;
        }

        if (yRight == 0.0 || yLeft * yRight < 0.0) {
            const double crossing = InterpolateCrossing(left, right, target);
            result.errHigh = crossing - result.points[result.bestIndex].fixedValue;
            break;
        }
    }

    if (result.fixedIndex >= 0 && result.fixedIndex < 3) {
        result.parErrLow[result.fixedIndex] = result.errLow;
        result.parErrHigh[result.fixedIndex] = result.errHigh;
    }
}

void FillBestPointMinosErrors(ScanResult &result,
                              TGraphErrors *gFit,
                              TGraphAsymmErrors *gr,
                              TH1F *eStat)
{
    result.minosOk = false;

    if (result.points.empty() || result.bestIndex < 0) {
        return;
    }

    const double fixedValue = result.points[result.bestIndex].fixedValue;

    TF1 *f = new TF1(("f_minos_" + result.name).c_str(), TsallisFuncAll, 0.0, 10.0, 3);
    ConfigureFunction(f);
    f->SetParameters(result.bestNor, result.bestQ, result.bestT);
    f->SetParameter(result.fixedIndex, fixedValue);

    for (int ipar = 0; ipar < 3; ++ipar) {
        if (ipar == result.fixedIndex) {
            f->FixParameter(ipar, fixedValue);
        } else {
            f->ReleaseParameter(ipar);
        }
    }

    TFitResultPtr r = gFit->Fit(f, "SQENR");
    if (int(r) != 0) {
        return;
    }

    result.bestNor = f->GetParameter(0);
    result.bestQ = f->GetParameter(1);
    result.bestT = f->GetParameter(2);
    result.bestChi2 = Chi2Stat(gr, eStat, f);

    result.points[result.bestIndex].Nor = result.bestNor;
    result.points[result.bestIndex].q = result.bestQ;
    result.points[result.bestIndex].T = result.bestT;
    result.points[result.bestIndex].chi2 = result.bestChi2;

    if (result.fixedIndex == 0) {
        result.bestNor = fixedValue;
    } else if (result.fixedIndex == 1) {
        result.bestQ = fixedValue;
    } else {
        result.bestT = fixedValue;
    }

    ComputeOneSigmaErrors(result);

    for (int ipar = 0; ipar < 3; ++ipar) {
        if (ipar == result.fixedIndex) {
            continue;
        }

        result.parErrLow[ipar] = std::fabs(r->LowerError(ipar));
        result.parErrHigh[ipar] = r->UpperError(ipar);
    }

    result.minosOk = true;
}

ScanResult RunProfileScan(const string &name,
                          const string &parameterName,
                          int fixedIndex,
                          double scanMin,
                          double scanMax,
                          int nScan,
                          TGraphErrors *gFit,
                          TGraphAsymmErrors *gr,
                          TH1F *eStat,
                          int nStatPoints)
{
    ScanResult result;
    result.name = name;
    result.parameterName = parameterName;
    result.fixedIndex = fixedIndex;
    result.bestIndex = -1;
    result.bestNor = 0.0;
    result.bestQ = 0.0;
    result.bestT = 0.0;
    result.bestChi2 = kBadChi2;
    result.errLow = numeric_limits<double>::quiet_NaN();
    result.errHigh = numeric_limits<double>::quiet_NaN();
    InitializeParameterErrors(result);
    result.minosOk = false;
    result.ndf = nStatPoints - 2;

    TF1 *f = new TF1(("f_" + name).c_str(), TsallisFuncAll, 0.0, 10.0, 3);
    ConfigureFunction(f);

    double initNor = 10.0;
    double initQ = 1.11;
    double initT = 0.10;

    for (int i = 0; i <= nScan; ++i) {
        const double fixedValue = scanMin + (scanMax - scanMin) * i / nScan;

        ScanPoint point;
        if (!FitAtFixedValue(gFit, gr, eStat, f, fixedIndex, fixedValue,
                             initNor, initQ, initT, point)) {
            continue;
        }

        result.points.push_back(point);

        if (point.chi2 < result.bestChi2) {
            result.bestIndex = (int)result.points.size() - 1;
            result.bestChi2 = point.chi2;
            result.bestNor = point.Nor;
            result.bestQ = point.q;
            result.bestT = point.T;
        }
    }

    ComputeOneSigmaErrors(result);
    FillBestPointMinosErrors(result, gFit, gr, eStat);

    return result;
}

string FormatDouble(double value, int precision = 8)
{
    if (!std::isfinite(value)) {
        return "NaN";
    }

    ostringstream out;
    out << setprecision(precision) << value;
    return out.str();
}

double FixedBestValue(const ScanResult &result)
{
    if (result.fixedIndex == 0) {
        return result.bestNor;
    }
    if (result.fixedIndex == 1) {
        return result.bestQ;
    }
    return result.bestT;
}

double ParameterValue(const ScanResult &result, int ipar)
{
    if (ipar == 0) {
        return result.bestNor;
    }
    if (ipar == 1) {
        return result.bestQ;
    }
    return result.bestT;
}

const char *ParameterName(int ipar)
{
    if (ipar == 0) {
        return "Nor";
    }
    if (ipar == 1) {
        return "q";
    }
    return "T";
}

void PrintParameterLine(ostream &out, const ScanResult &result, int ipar)
{
    const char *source = (ipar == result.fixedIndex) ? "profile" : "MINOS";

    out << left << setw(3) << ParameterName(ipar) << right
        << " = " << FormatDouble(ParameterValue(result, ipar))
        << " +" << FormatDouble(result.parErrHigh[ipar])
        << " -" << FormatDouble(result.parErrLow[ipar])
        << "  (" << source << ")\n";
}

void PrintScanResult(ostream &out, const ScanResult &result)
{
    out << "\n=== " << result.name << " ===\n";
    PrintParameterLine(out, result, 0);
    PrintParameterLine(out, result, 1);
    PrintParameterLine(out, result, 2);
    out << "chi2_stat = " << FormatDouble(result.bestChi2) << "\n";
    out << "ndf = " << result.ndf << "\n";
    out << "chi2_stat/ndf = " << FormatDouble(result.bestChi2 / result.ndf) << "\n";

    if (!std::isfinite(result.errLow) || !std::isfinite(result.errHigh)) {
        out << "warning: chi2_min + 1 crossing was not found on at least one side.\n";
    }

    if (!result.minosOk) {
        out << "warning: MINOS failed for the two non-fixed parameters.\n";
    }
}

void StyleDataGraph(TGraphAsymmErrors *gr)
{
    gr->SetLineColor(kCyan - 2);
    gr->SetLineWidth(2);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(1);
    gr->GetXaxis()->SetLimits(0, 12);
    gr->SetMinimum(TMath::Power(10, -8));
    gr->SetMaximum(TMath::Power(10, 2));
    gr->SetTitle("");
    gr->GetXaxis()->CenterTitle();
    gr->GetXaxis()->SetTitle("p_{T} (GeV)");
    gr->GetXaxis()->SetTitleOffset(1.3);
    gr->GetYaxis()->CenterTitle();
    gr->GetYaxis()->SetLabelSize(0.03);
    gr->GetYaxis()->SetTitle("1/N 1/(2#pip_{T}) d^{2}N/(dydp_{T}) (GeV/c)^{-2}");
    gr->GetYaxis()->SetTitleOffset(1.3);
}

void SaveProfilePlot(const ScanResult &result, const char *fileName)
{
    if (result.points.empty()) {
        cerr << "No valid scan points for " << result.name << endl;
        return;
    }

    TCanvas *c = new TCanvas(("c_" + result.name).c_str(), result.name.c_str(), 800, 600);
    c->SetTickx();
    c->SetTicky();

    TGraph *g = new TGraph();
    g->SetName(("g_" + result.name).c_str());

    double yMin = kBadChi2;
    double yMax = -kBadChi2;
    for (int i = 0; i < (int)result.points.size(); ++i) {
        g->SetPoint(i, result.points[i].fixedValue, result.points[i].chi2);
        yMin = std::min(yMin, result.points[i].chi2);
        yMax = std::max(yMax, result.points[i].chi2);
    }

    const double target = result.bestChi2 + 1.0;
    yMax = std::min(yMax, result.bestChi2 + 50.0);
    if (yMax <= target + 0.5) {
        yMax = target + 1.0;
    }

    g->SetTitle(("#chi^{2}_{stat} vs " + result.parameterName + ";" +
                 result.parameterName + ";#chi^{2}_{stat}")
                    .c_str());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.6);
    g->SetLineWidth(2);
    g->SetMinimum(yMin - 0.05 * (yMax - yMin));
    g->SetMaximum(yMax);
    g->Draw("ALP");

    const double bestValue = FixedBestValue(result);
    TLine *lineBest = new TLine(bestValue, yMin, bestValue, yMax);
    lineBest->SetLineColor(kRed + 1);
    lineBest->SetLineStyle(2);
    lineBest->SetLineWidth(2);
    lineBest->Draw("same");

    TLine *lineSigma = new TLine(result.points.front().fixedValue, target,
                                 result.points.back().fixedValue, target);
    lineSigma->SetLineColor(kBlue + 1);
    lineSigma->SetLineStyle(7);
    lineSigma->SetLineWidth(2);
    lineSigma->Draw("same");

    if (std::isfinite(result.errLow)) {
        const double xLow = bestValue - result.errLow;
        TLine *lineLow = new TLine(xLow, yMin, xLow, target);
        lineLow->SetLineColor(kGreen + 2);
        lineLow->SetLineStyle(3);
        lineLow->Draw("same");
    }

    if (std::isfinite(result.errHigh)) {
        const double xHigh = bestValue + result.errHigh;
        TLine *lineHigh = new TLine(xHigh, yMin, xHigh, target);
        lineHigh->SetLineColor(kGreen + 2);
        lineHigh->SetLineStyle(3);
        lineHigh->Draw("same");
    }

    TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g, "#chi^{2}_{stat}", "lp");
    leg->AddEntry(lineBest, "best", "l");
    leg->AddEntry(lineSigma, "#chi^{2}_{min}+1", "l");
    leg->Draw("same");

    c->SaveAs(OutputPath(fileName).c_str());
}

void SaveBestFitPlot(TGraphAsymmErrors *gr,
                     double Nor,
                     double q,
                     double T,
                     const string &title,
                     const char *fileName)
{
    TCanvas *c = new TCanvas(("c_fit_" + title).c_str(), title.c_str(), 550, 500);
    c->SetTickx();
    c->SetTicky();
    c->SetLogy();

    StyleDataGraph(gr);
    gr->Draw("AP");

    TF1 *f = new TF1(("f_fit_" + title).c_str(), TsallisFuncAll, 0.0, 10.0, 3);
    ConfigureFunction(f);
    f->SetParameters(Nor, q, T);
    f->SetLineColor(kRed + 1);
    f->SetLineWidth(2);
    f->Draw("same");

    TLegend *leg = new TLegend(0.42, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr, "HEPData", "p");
    leg->AddEntry(f, "Tsallis fit", "l");
    leg->AddEntry((TObject*)0, Form("Nor = %.5g", Nor), "");
    leg->AddEntry((TObject*)0, Form("q = %.6f, T = %.6f", q, T), "");
    leg->Draw("same");

    c->SaveAs(OutputPath(fileName).c_str());
}

int FreeParameterIndex(int xIndex, int yIndex)
{
    for (int ipar = 0; ipar < 3; ++ipar) {
        if (ipar != xIndex && ipar != yIndex) {
            return ipar;
        }
    }

    return -1;
}

double DefaultParameterValue(int ipar)
{
    if (ipar == 0) {
        return 10.0;
    }
    if (ipar == 1) {
        return 1.11;
    }
    return 0.10;
}

PairProfileResult SavePairProfilePlots(TGraphErrors *gFit,
                                       TGraphAsymmErrors *gr,
                                       TH1F *eStat,
                                       int nStatPoints,
                                       const string &name,
                                       int xIndex,
                                       int yIndex,
                                       double xMin,
                                       double xMax,
                                       int nX,
                                       double yMin,
                                       double yMax,
                                       int nY)
{
    PairProfileResult result;
    result.name = name;
    result.xName = ParameterName(xIndex);
    result.yName = ParameterName(yIndex);
    result.xIndex = xIndex;
    result.yIndex = yIndex;
    result.freeIndex = FreeParameterIndex(xIndex, yIndex);
    result.bestChi2 = kBadChi2;
    result.ndf = nStatPoints - 1;

    for (int ipar = 0; ipar < 3; ++ipar) {
        result.bestPars[ipar] = 0.0;
    }
    InitializePairProfileErrors(result);

    TF1 *f = new TF1(("f_pair_" + name).c_str(), TsallisFuncAll, 0.0, 10.0, 3);
    ConfigureFunction(f);

    TH2D *hChi2 = new TH2D(("hChi2_" + name).c_str(),
                           Form("#chi^{2}_{stat}(%s,%s);%s;%s",
                                result.xName.c_str(), result.yName.c_str(),
                                result.xName.c_str(), result.yName.c_str()),
                           nX, xMin, xMax, nY, yMin, yMax);

    vector<vector<PairProfilePoint> > profileGrid(nX, vector<PairProfilePoint>(nY));

    double initPars[3] = {DefaultParameterValue(0),
                          DefaultParameterValue(1),
                          DefaultParameterValue(2)};

    for (int ix = 0; ix < nX; ++ix) {
        const double xValue = xMin + (xMax - xMin) * ix / (nX - 1);

        for (int iy = 0; iy < nY; ++iy) {
            const double yValue = yMin + (yMax - yMin) * iy / (nY - 1);

            double pars[3] = {initPars[0], initPars[1], initPars[2]};
            pars[xIndex] = xValue;
            pars[yIndex] = yValue;

            f->SetParameters(pars);
            ConfigureFunction(f);

            for (int ipar = 0; ipar < 3; ++ipar) {
                if (ipar == xIndex) {
                    f->FixParameter(ipar, xValue);
                } else if (ipar == yIndex) {
                    f->FixParameter(ipar, yValue);
                } else {
                    f->ReleaseParameter(ipar);
                }
            }

            TFitResultPtr fitResult = gFit->Fit(f, "SQNR");
            if (int(fitResult) != 0 ||
                f->GetParameter(0) <= 0 ||
                f->GetParameter(1) <= 1.0 ||
                f->GetParameter(2) <= 0) {
                hChi2->SetBinContent(ix + 1, iy + 1, kBadChi2);
                continue;
            }

            const double chi2 = Chi2Stat(gr, eStat, f);
            hChi2->SetBinContent(ix + 1, iy + 1, chi2);

            PairProfilePoint &profilePoint = profileGrid[ix][iy];
            profilePoint.valid = true;
            profilePoint.chi2 = chi2;
            for (int ipar = 0; ipar < 3; ++ipar) {
                profilePoint.pars[ipar] = f->GetParameter(ipar);
            }
            profilePoint.pars[xIndex] = xValue;
            profilePoint.pars[yIndex] = yValue;

            if (chi2 < result.bestChi2) {
                result.bestChi2 = chi2;
                result.bestPars[0] = f->GetParameter(0);
                result.bestPars[1] = f->GetParameter(1);
                result.bestPars[2] = f->GetParameter(2);
                result.bestPars[xIndex] = xValue;
                result.bestPars[yIndex] = yValue;
            }

            if (result.freeIndex >= 0) {
                initPars[result.freeIndex] = f->GetParameter(result.freeIndex);
            }
        }
    }

    ComputePairProfileErrors(result, profileGrid);

    gStyle->SetOptStat(0);

    TH2D *hDelta = (TH2D*)hChi2->Clone(("hDelta_" + name).c_str());
    hDelta->SetTitle(Form("#Delta#chi^{2}_{stat}(%s,%s);%s;%s",
                          result.xName.c_str(), result.yName.c_str(),
                          result.xName.c_str(), result.yName.c_str()));

    for (int ix = 1; ix <= hDelta->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hDelta->GetNbinsY(); ++iy) {
            const double v = hDelta->GetBinContent(ix, iy);
            if (v <= 0 || v > result.bestChi2 + 500.0) {
                hDelta->SetBinContent(ix, iy, 500.0);
            } else {
                hDelta->SetBinContent(ix, iy, v - result.bestChi2);
            }
        }
    }

    double deltaLevels[6] = {1.0, 2.30, 4.0, 6.18, 9.0, 11.83};

    TCanvas *cEllipse = new TCanvas(("c_ellipse_" + name).c_str(),
                                    ("chi2 ellipse " + name).c_str(), 650, 550);
    hDelta->SetMinimum(0.0);
    hDelta->SetMaximum(20.0);
    hDelta->Draw("COLZ");

    TH2D *hEllipseLines = (TH2D*)hDelta->Clone(("hEllipseLines_" + name).c_str());
    hEllipseLines->SetContour(6, deltaLevels);
    hEllipseLines->SetLineColor(kBlack);
    hEllipseLines->SetLineWidth(2);
    hEllipseLines->Draw("CONT3 SAME");

    TMarker *mBest = new TMarker(result.bestPars[xIndex], result.bestPars[yIndex], 29);
    mBest->SetMarkerColor(kRed);
    mBest->SetMarkerSize(2.0);
    mBest->Draw("same");
    cEllipse->SaveAs(OutputPath("chi2_" + name + "_ellipse.png").c_str());

    TCanvas *cContour = new TCanvas(("c_contour_" + name).c_str(),
                                    ("chi2 contour " + name).c_str(), 650, 550);
    double absoluteLevels[6] = {
        result.bestChi2 + 1.0,
        result.bestChi2 + 2.30,
        result.bestChi2 + 4.0,
        result.bestChi2 + 6.18,
        result.bestChi2 + 9.0,
        result.bestChi2 + 11.83
    };

    hChi2->SetContour(6, absoluteLevels);
    hChi2->Draw("CONT1");

    TMarker *mBest2 = new TMarker(result.bestPars[xIndex], result.bestPars[yIndex], 29);
    mBest2->SetMarkerColor(kRed);
    mBest2->SetMarkerSize(2.0);
    mBest2->Draw("same");
    cContour->SaveAs(OutputPath("chi2_" + name + "_contour.png").c_str());

    return result;
}

ContourResult SaveQTContourPlots(TGraphErrors *gFit,
                                 TGraphAsymmErrors *gr,
                                 TH1F *eStat,
                                 int nStatPoints)
{
    ContourResult result;
    result.bestNor = 0.0;
    result.bestQ = 0.0;
    result.bestT = 0.0;
    result.bestChi2 = kBadChi2;
    result.ndf = nStatPoints - 1;

    TF1 *f = new TF1("f_qT_contour_all", TsallisFuncAll, 0.0, 10.0, 3);
    ConfigureFunction(f);

    TH2D *hChi2 = new TH2D("hChi2_qT_all",
                           "#chi^{2}_{stat}(q,T);q;T",
                           kNQContour, kQMin, kQMax,
                           kNTContour, kTMin, kTMax);

    double initNor = 10.0;

    for (int iq = 0; iq < kNQContour; ++iq) {
        const double qValue = kQMin + (kQMax - kQMin) * iq / (kNQContour - 1);

        for (int iT = 0; iT < kNTContour; ++iT) {
            const double TValue = kTMin + (kTMax - kTMin) * iT / (kNTContour - 1);

            f->SetParameters(initNor, qValue, TValue);
            ConfigureFunction(f);
            f->ReleaseParameter(0);
            f->FixParameter(1, qValue);
            f->FixParameter(2, TValue);

            TFitResultPtr fitResult = gFit->Fit(f, "SQNR");
            if (int(fitResult) != 0 || f->GetParameter(0) <= 0) {
                hChi2->SetBinContent(iq + 1, iT + 1, kBadChi2);
                continue;
            }

            const double Nor = f->GetParameter(0);
            const double chi2 = Chi2Stat(gr, eStat, f);
            hChi2->SetBinContent(iq + 1, iT + 1, chi2);

            if (chi2 < result.bestChi2) {
                result.bestChi2 = chi2;
                result.bestNor = Nor;
                result.bestQ = qValue;
                result.bestT = TValue;
            }

            initNor = Nor;
        }
    }

    gStyle->SetOptStat(0);

    TCanvas *cColz = new TCanvas("c_qT_colz_all", "chi2 q T colz", 550, 500);
    TH2D *hDelta = (TH2D*)hChi2->Clone("hDelta_qT_all");
    hDelta->SetTitle("#Delta#chi^{2}_{stat}(q,T);q;T");

    for (int ix = 1; ix <= hDelta->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hDelta->GetNbinsY(); ++iy) {
            const double v = hDelta->GetBinContent(ix, iy);
            if (v <= 0 || v > result.bestChi2 + 500.0) {
                hDelta->SetBinContent(ix, iy, 500.0);
            } else {
                hDelta->SetBinContent(ix, iy, v - result.bestChi2);
            }
        }
    }

    hDelta->SetMinimum(0.0);
    hDelta->SetMaximum(20.0);
    hDelta->Draw("COLZ");

    TMarker *mBest = new TMarker(result.bestQ, result.bestT, 29);
    mBest->SetMarkerColor(kRed);
    mBest->SetMarkerSize(2.0);
    mBest->Draw("same");
    cColz->SaveAs(OutputPath("chi2_qT_colz.png").c_str());

    TCanvas *cContour = new TCanvas("c_qT_contour_all", "chi2 q T contour", 550, 500);
    double levels[6] = {
        result.bestChi2 + 1.0,
        result.bestChi2 + 2.30,
        result.bestChi2 + 4.0,
        result.bestChi2 + 6.18,
        result.bestChi2 + 9.0,
        result.bestChi2 + 11.83
    };

    hChi2->SetContour(6, levels);
    hChi2->Draw("CONT1");

    TMarker *mBest2 = new TMarker(result.bestQ, result.bestT, 29);
    mBest2->SetMarkerColor(kRed);
    mBest2->SetMarkerSize(2.0);
    mBest2->Draw("same");
    cContour->SaveAs(OutputPath("chi2_qT_contour.png").c_str());

    SaveBestFitPlot(gr, result.bestNor, result.bestQ, result.bestT,
                    "qT_contour", "best_fit_qT.png");

    return result;
}

void PrintContourResult(ostream &out, const ContourResult &result)
{
    out << "\n=== q-T contour, fixed q and T, fit Nor ===\n";
    out << "Nor = " << FormatDouble(result.bestNor) << "\n";
    out << "q   = " << FormatDouble(result.bestQ) << "\n";
    out << "T   = " << FormatDouble(result.bestT) << "\n";
    out << "chi2_stat = " << FormatDouble(result.bestChi2) << "\n";
    out << "ndf = " << result.ndf << "\n";
    out << "chi2_stat/ndf = " << FormatDouble(result.bestChi2 / result.ndf) << "\n";
}

void PrintPairParameterLine(ostream &out, const PairProfileResult &result, int ipar)
{
    out << left << setw(3) << ParameterName(ipar) << right
        << " = " << FormatDouble(result.bestPars[ipar])
        << " +" << FormatDouble(result.parErrHigh[ipar])
        << " -" << FormatDouble(result.parErrLow[ipar])
        << "\n";
}

void PrintPairProfileResult(ostream &out, const PairProfileResult &result)
{
    out << "\n=== " << result.name << " 2D profile contour ===\n";
    out << "fixed axes: " << result.xName << ", " << result.yName
        << "; fitted parameter: " << ParameterName(result.freeIndex) << "\n";
    out << "errors: Delta chi2_stat = " << FormatDouble(result.errorDeltaChi2)
        << " profile-contour projection\n";
    PrintPairParameterLine(out, result, 0);
    PrintPairParameterLine(out, result, 1);
    PrintPairParameterLine(out, result, 2);
    out << "chi2_stat = " << FormatDouble(result.bestChi2) << "\n";
    out << "ndf = " << result.ndf << "\n";
    out << "chi2_stat/ndf = " << FormatDouble(result.bestChi2 / result.ndf) << "\n";

    if (!std::isfinite(result.parErrLow[0]) ||
        !std::isfinite(result.parErrHigh[0]) ||
        !std::isfinite(result.parErrLow[1]) ||
        !std::isfinite(result.parErrHigh[1]) ||
        !std::isfinite(result.parErrLow[2]) ||
        !std::isfinite(result.parErrHigh[2])) {
        out << "warning: Delta chi2 crossing was not found on at least one side.\n";
    }
}

void tsallis_graph1d_y2_scans()
{
    gErrorIgnoreLevel = kFatal;
    gSystem->mkdir(kOutputDir.c_str(), true);

    TFile *file = TFile::Open("HEPData-ins2061074-v1-Figure_15.root");
    if (!file || file->IsZombie()) {
        cerr << "Failed to open HEPData-ins2061074-v1-Figure_15.root" << endl;
        return;
    }

    TDirectoryFile *dir = (TDirectoryFile*)file->Get("Figure 15");
    if (!dir) {
        cerr << "Failed to find directory: Figure 15" << endl;
        return;
    }

    TGraphAsymmErrors *gr = (TGraphAsymmErrors*)dir->Get("Graph1D_y2");
    TH1F *eStat = (TH1F*)dir->Get("Hist1D_y2_e1");
    TH1F *eSys = (TH1F*)dir->Get("Hist1D_y2_e2");

    if (!gr || !eStat || !eSys) {
        cerr << "Failed to load Graph1D_y2, Hist1D_y2_e1, or Hist1D_y2_e2" << endl;
        return;
    }

    const int nStatPoints = CountStatPoints(gr, eStat);
    TGraphErrors *gFit = BuildStatFitGraph(gr, eStat);

    ScanResult fixedNor = RunProfileScan("fixed_Nor", "Nor", 0,
                                         kNorMin, kNorMax, kScan1D,
                                         gFit, gr, eStat, nStatPoints);
    ScanResult fixedQ = RunProfileScan("fixed_q", "q", 1,
                                       kQMin, kQMax, kScan1D,
                                       gFit, gr, eStat, nStatPoints);
    ScanResult fixedT = RunProfileScan("fixed_T", "T", 2,
                                       kTMin, kTMax, kScan1D,
                                       gFit, gr, eStat, nStatPoints);

    SaveProfilePlot(fixedNor, "chi2_vs_Nor.png");
    SaveProfilePlot(fixedQ, "chi2_vs_q.png");
    SaveProfilePlot(fixedT, "chi2_vs_T.png");

    SaveBestFitPlot(gr, fixedNor.bestNor, fixedNor.bestQ, fixedNor.bestT,
                    "fixed_Nor", "best_fit_fixed_Nor.png");
    SaveBestFitPlot(gr, fixedQ.bestNor, fixedQ.bestQ, fixedQ.bestT,
                    "fixed_q", "best_fit_fixed_q.png");
    SaveBestFitPlot(gr, fixedT.bestNor, fixedT.bestQ, fixedT.bestT,
                    "fixed_T", "best_fit_fixed_T.png");

    PairProfileResult qTContour = SavePairProfilePlots(gFit, gr, eStat, nStatPoints,
                                                       "qT", 1, 2,
                                                       kQMin, kQMax, kNQContour,
                                                       kTMin, kTMax, kNTContour);
    PairProfileResult NorQContour = SavePairProfilePlots(gFit, gr, eStat, nStatPoints,
                                                         "Nor_q", 0, 1,
                                                         kNorMin, kNorMax, kNQContour,
                                                         kQMin, kQMax, kNQContour);
    PairProfileResult NorTContour = SavePairProfilePlots(gFit, gr, eStat, nStatPoints,
                                                         "Nor_T", 0, 2,
                                                         kNorMin, kNorMax, kNTContour,
                                                         kTMin, kTMax, kNTContour);

    SaveBestFitPlot(gr, qTContour.bestPars[0], qTContour.bestPars[1], qTContour.bestPars[2],
                    "qT_contour", "best_fit_qT.png");
    SaveBestFitPlot(gr, NorQContour.bestPars[0], NorQContour.bestPars[1], NorQContour.bestPars[2],
                    "Nor_q_contour", "best_fit_Nor_q.png");
    SaveBestFitPlot(gr, NorTContour.bestPars[0], NorTContour.bestPars[1], NorTContour.bestPars[2],
                    "Nor_T_contour", "best_fit_Nor_T.png");

    ofstream resultFile(OutputPath("fit_results.txt").c_str());

    cout << setprecision(10);
    resultFile << setprecision(10);

    cout << "Input: Figure 15 / Graph1D_y2" << endl;
    cout << "chi2: statistical errors from Hist1D_y2_e1" << endl;
    cout << "N(stat points) = " << nStatPoints << endl;

    resultFile << "Input: Figure 15 / Graph1D_y2\n";
    resultFile << "chi2: statistical errors from Hist1D_y2_e1\n";
    resultFile << "N(stat points) = " << nStatPoints << "\n";

    PrintScanResult(cout, fixedNor);
    PrintScanResult(cout, fixedQ);
    PrintScanResult(cout, fixedT);
    PrintPairProfileResult(cout, qTContour);
    PrintPairProfileResult(cout, NorQContour);
    PrintPairProfileResult(cout, NorTContour);

    PrintScanResult(resultFile, fixedNor);
    PrintScanResult(resultFile, fixedQ);
    PrintScanResult(resultFile, fixedT);
    PrintPairProfileResult(resultFile, qTContour);
    PrintPairProfileResult(resultFile, NorQContour);
    PrintPairProfileResult(resultFile, NorTContour);

    cout << "\nSaved output directory: " << kOutputDir << endl;
    cout << "Saved plots:" << endl;
    cout << "  " << OutputPath("chi2_vs_Nor.png") << endl;
    cout << "  " << OutputPath("chi2_vs_q.png") << endl;
    cout << "  " << OutputPath("chi2_vs_T.png") << endl;
    cout << "  " << OutputPath("chi2_qT_ellipse.png") << endl;
    cout << "  " << OutputPath("chi2_qT_contour.png") << endl;
    cout << "  " << OutputPath("chi2_Nor_q_ellipse.png") << endl;
    cout << "  " << OutputPath("chi2_Nor_q_contour.png") << endl;
    cout << "  " << OutputPath("chi2_Nor_T_ellipse.png") << endl;
    cout << "  " << OutputPath("chi2_Nor_T_contour.png") << endl;
    cout << "  " << OutputPath("best_fit_fixed_Nor.png") << endl;
    cout << "  " << OutputPath("best_fit_fixed_q.png") << endl;
    cout << "  " << OutputPath("best_fit_fixed_T.png") << endl;
    cout << "  " << OutputPath("best_fit_qT.png") << endl;
    cout << "  " << OutputPath("best_fit_Nor_q.png") << endl;
    cout << "  " << OutputPath("best_fit_Nor_T.png") << endl;
    cout << "Saved parameters: " << OutputPath("fit_results.txt") << endl;
}
