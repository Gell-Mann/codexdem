// Tsallis blast-wave fit for the pion pT spectrum in HEPData Figure 15.
// Run with: root -l -b -q TBW.C

#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TStyle.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {

constexpr double kPionMassGeV = 0.13957039;
constexpr double kRadius = 22.0;
constexpr int kRadialSteps = 80;
constexpr int kPhiSteps = 120;

double TsallisBlastWave(double *x, double *par)
{
    const double pT = x[0];
    const double T0 = par[0];
    const double n0 = par[1];
    const double betaS = par[2];
    const double q = par[3];
    const double norm = par[4];

    if (T0 <= 0.0 || q <= 1.0 || betaS < 0.0 || betaS >= 1.0 || norm < 0.0) {
        return 0.0;
    }

    const double mt = TMath::Sqrt(pT * pT + kPionMassGeV * kPionMassGeV);
    const double dr = kRadius / kRadialSteps;
    const double dphi = 2.0 * TMath::Pi() / kPhiSteps;
    double integral = 0.0;

    // Deterministic midpoint integration avoids MC noise inside the minimizer.
    for (int ir = 0; ir < kRadialSteps; ++ir) {
        const double r = (ir + 0.5) * dr;
        const double betaR = betaS * TMath::Power(r / kRadius, n0);
        if (betaR >= 1.0) {
            continue;
        }

        const double rho = TMath::ATanH(betaR);
        const double sinhRho = TMath::SinH(rho);
        const double coshRho = TMath::CosH(rho);

        for (int iphi = 0; iphi < kPhiSteps; ++iphi) {
            const double phi = -TMath::Pi() + (iphi + 0.5) * dphi;
            const double energyTerm = mt * coshRho - pT * sinhRho * TMath::Cos(phi);
            const double base = 1.0 + (q - 1.0) * energyTerm / T0;
            if (base > 0.0) {
                integral += r * TMath::Power(base, -1.0 / (q - 1.0));
            }
        }
    }

    integral *= dr * dphi;
    return norm * mt * integral;
}

void SaveFitSummary(const TF1 &fit, const TFitResultPtr &result, const char *path)
{
    std::ofstream out(path);
    out << std::setprecision(10);
    out << "name,value,error,lower_limit,upper_limit\n";
    for (int ipar = 0; ipar < fit.GetNpar(); ++ipar) {
        double low = 0.0;
        double high = 0.0;
        fit.GetParLimits(ipar, low, high);
        out << fit.GetParName(ipar) << ','
            << fit.GetParameter(ipar) << ','
            << fit.GetParError(ipar) << ','
            << low << ','
            << high << '\n';
    }
    out << "chi2," << fit.GetChisquare() << ",,,\n";
    out << "ndf," << fit.GetNDF() << ",,,\n";
    out << "chi2_ndf," << (fit.GetNDF() > 0 ? fit.GetChisquare() / fit.GetNDF() : 0.0) << ",,,\n";
    out << "fit_status," << static_cast<int>(result) << ",,,\n";
}

void SaveResiduals(const TGraphAsymmErrors &graph, const TF1 &fit, const char *path)
{
    std::ofstream out(path);
    out << std::setprecision(10);
    out << "pT,y,fit,y_error_avg,residual,pull\n";
    for (int i = 0; i < graph.GetN(); ++i) {
        double pT = 0.0;
        double y = 0.0;
        graph.GetPoint(i, pT, y);
        const double yFit = fit.Eval(pT);
        const double yErr = 0.5 * (graph.GetErrorYlow(i) + graph.GetErrorYhigh(i));
        const double residual = y - yFit;
        const double pull = yErr > 0.0 ? residual / yErr : 0.0;
        out << pT << ','
            << y << ','
            << yFit << ','
            << yErr << ','
            << residual << ','
            << pull << '\n';
    }
}

} // namespace

void TBW()
{
    constexpr const char *inputPath = "HEPData-ins2061074-v1-Figure_15.root";
    constexpr const char *directoryName = "Figure 15";
    constexpr const char *graphName = "Graph1D_y1";
    constexpr double xMin = 0.0;
    constexpr double xMax = 12.0;

    const int oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
    if (!input || input->IsZombie()) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot open " << inputPath << '\n';
        return;
    }

    auto *directory = dynamic_cast<TDirectoryFile *>(input->Get(directoryName));
    if (!directory) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot find directory \"" << directoryName << "\"\n";
        return;
    }

    auto *inputGraph = dynamic_cast<TGraphAsymmErrors *>(directory->Get(graphName));
    if (!inputGraph) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot find graph \"" << graphName << "\"\n";
        return;
    }
    std::unique_ptr<TGraphAsymmErrors> graph(
        dynamic_cast<TGraphAsymmErrors *>(inputGraph->Clone("pionSpectrum")));
    if (!graph) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot clone graph \"" << graphName << "\"\n";
        return;
    }
    input.reset();
    gErrorIgnoreLevel = oldErrorIgnoreLevel;

    TF1 fit("tbwFit", TsallisBlastWave, xMin, xMax, 5);
    fit.SetLineColor(kRed + 1);
    fit.SetLineWidth(2);
    fit.SetParNames("T0", "n0", "betaS", "q", "norm");
    fit.SetParameters(0.10, 2.0, 0.10, 1.10, 1.0e-3);
    fit.SetParLimits(0, 0.05, 0.20);
    fit.FixParameter(1, 2.0);
    fit.FixParameter(2, 0.10);
    fit.SetParLimits(3, 1.0001, 1.20);
    fit.SetParLimits(4, 0.0, 1.0e6);

    TCanvas canvas("fit1", "TBW fit", 900, 650);
    canvas.SetTickx();
    canvas.SetTicky();
    canvas.SetLogy();

    graph->SetLineColor(kCyan - 2);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue + 1);
    graph->SetMarkerSize(1);
    graph->GetXaxis()->SetLimits(xMin, xMax);
    graph->SetMinimum(1.0e-9);
    graph->SetMaximum(1.0e3);

    gStyle->SetOptFit(1111);
    gStyle->SetTitleY(0.96);
    graph->SetTitle("TBW fit #pi^{+}");
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetLabelSize(0.03);
    graph->GetYaxis()->SetTitle("(1/2#pip_{T}) d^{2}N/dydp_{T} [(GeV/c)^{-2}]");
    graph->GetYaxis()->SetTitleOffset(1.35);

    graph->Draw("AP");
    gErrorIgnoreLevel = kFatal;
    TFitResultPtr result = graph->Fit(&fit, "RS0Q");
    fit.Draw("same");

    canvas.SaveAs("TBW_fit.pdf");
    canvas.SaveAs("TBW_fit.png");
    SaveFitSummary(fit, result, "TBW_fit_summary.csv");
    SaveResiduals(*graph, fit, "TBW_fit_residuals.csv");

    std::cout << "Fit status = " << static_cast<int>(result) << '\n';
    for (int ipar = 0; ipar < fit.GetNpar(); ++ipar) {
        double low = 0.0;
        double high = 0.0;
        fit.GetParLimits(ipar, low, high);
        std::cout << fit.GetParName(ipar) << " = " << fit.GetParameter(ipar)
                  << " +/- " << fit.GetParError(ipar);
        if (low == high && low == fit.GetParameter(ipar)) {
            std::cout << " (fixed)";
        }
        std::cout << '\n';
    }
    std::cout << "chi2/ndf = " << fit.GetChisquare() << '/' << fit.GetNDF();
    if (fit.GetNDF() > 0) {
        std::cout << " = " << fit.GetChisquare() / fit.GetNDF();
    }
    std::cout << '\n';
    std::cout << "Saved: TBW_fit.pdf, TBW_fit.png, TBW_fit_summary.csv, TBW_fit_residuals.csv\n";
}
