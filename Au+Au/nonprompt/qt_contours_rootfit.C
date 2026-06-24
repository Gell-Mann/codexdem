#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TError.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

struct RootFitDataSet {
    const char *graphName;
    const char *statPlusName;
    const char *statMinusName;
    const char *sysPlusName;
    const char *sysMinusName;
    const char *label;
    int color;
    std::shared_ptr<TGraphErrors> graph;
};

struct RootFitResult {
    double nor = 0.;
    double q = 0.;
    double temp = 0.;
    double norError = 0.;
    double qError = 0.;
    double tempError = 0.;
    double chi2 = 0.;
    int ndf = 0;
    int status = -1;
    int covarianceStatus = -1;
};

double rootfit_shape(double *x, double *par)
{
    // par[0] = Nor, par[1] = q, par[2] = T。
    const double b = x[0];
    const double q = par[1];
    const double temp = par[2];
    const double z1 = TMath::Power(1. + (q - 1.) * b / temp,
                                   q / (1. - q));
    const double z2 = 2. * b / TMath::Power(2. * TMath::Pi(), 3);
    return par[0] * z1 * z2 * TMath::Power(10., 4);
}

double rootfit_symmetric_error(const TH1F *plus, const TH1F *minus,
                               int bin)
{
    return 0.5 * (std::abs(plus->GetBinContent(bin)) +
                  std::abs(minus->GetBinContent(bin)));
}

bool rootfit_load_dataset(TDirectoryFile *directory,
                          RootFitDataSet &dataSet,
                          double xMin, double xMax)
{
    auto *source =
        dynamic_cast<TGraphAsymmErrors *>(directory->Get(dataSet.graphName));
    auto *statPlus =
        dynamic_cast<TH1F *>(directory->Get(dataSet.statPlusName));
    auto *statMinus =
        dynamic_cast<TH1F *>(directory->Get(dataSet.statMinusName));
    auto *sysPlus =
        dynamic_cast<TH1F *>(directory->Get(dataSet.sysPlusName));
    auto *sysMinus =
        dynamic_cast<TH1F *>(directory->Get(dataSet.sysMinusName));
    if (!source || !statPlus || !statMinus || !sysPlus || !sysMinus) {
        std::cerr << "Error: cannot load objects for " << dataSet.label << '\n';
        return false;
    }

    dataSet.graph = std::make_shared<TGraphErrors>();
    dataSet.graph->SetName(Form("rootfit_data_%s", dataSet.label));
    for (int i = 0; i < source->GetN(); ++i) {
        const double x = source->GetPointX(i);
        if (x < xMin || x > xMax) continue;

        const int bin = i + 1;
        const double stat =
            rootfit_symmetric_error(statPlus, statMinus, bin);
        const double sys =
            rootfit_symmetric_error(sysPlus, sysMinus, bin);
        const double totalError = std::sqrt(stat * stat + sys * sys);
        if (totalError <= 0.) {
            std::cerr << "Warning: skip point with zero uncertainty: "
                      << dataSet.label << ", x=" << x << '\n';
            continue;
        }

        const int point = dataSet.graph->GetN();
        dataSet.graph->SetPoint(point, x, source->GetPointY(i));
        dataSet.graph->SetPointError(point, 0., totalError);
    }
    return dataSet.graph->GetN() > 0;
}

bool rootfit_dataset(RootFitDataSet &dataSet, double xMin, double xMax,
                     double norMax, RootFitResult &output,
                     std::unique_ptr<TGraph> &contour)
{
    auto function = std::make_unique<TF1>(
        Form("rootfit_function_%s", dataSet.label),
        rootfit_shape, xMin, xMax, 3);
    function->SetParNames("Nor", "q", "T");
    function->SetParameters(1., 1.09, 0.12);
    function->SetParLimits(0, 0., norMax);
    function->SetParLimits(1, 1.001, 1.5);
    function->SetParLimits(2, 0.001, 0.30);

    // S: 保存拟合及参数误差；R: 使用 TF1 范围；Q0: 静默且不自动绘图。
    // χ²和 Hessian/协方差拟合误差均由 TGraphErrors::Fit/Minuit 计算。
    TFitResultPtr fitResult = dataSet.graph->Fit(function.get(), "SRQ0");
    if (!fitResult.Get()) {
        std::cerr << "Error: ROOT did not return a fit result for "
                  << dataSet.label << '\n';
        return false;
    }

    output.status = fitResult->Status();
    output.covarianceStatus = fitResult->CovMatrixStatus();
    output.nor = fitResult->Parameter(0);
    output.q = fitResult->Parameter(1);
    output.temp = fitResult->Parameter(2);
    output.norError = fitResult->ParError(0);
    output.qError = fitResult->ParError(1);
    output.tempError = fitResult->ParError(2);
    output.chi2 = fitResult->Chi2();
    output.ndf = fitResult->Ndf();

    // 使用 ROOT 拟合返回的 q-T 协方差构造 68.3% 拟合误差椭圆。
    const TMatrixDSym covariance = fitResult->GetCovarianceMatrix();
    const double cqq = covariance(1, 1);
    const double cqt = covariance(1, 2);
    const double ctt = covariance(2, 2);
    if (cqq <= 0.) {
        std::cerr << "Error: invalid q variance for " << dataSet.label << '\n';
        return false;
    }
    const double l11 = std::sqrt(cqq);
    const double l21 = cqt / l11;
    const double l22Squared = ctt - l21 * l21;
    if (l22Squared <= 0.) {
        std::cerr << "Error: q-T covariance is not positive definite for "
                  << dataSet.label << '\n';
        return false;
    }
    const double l22 = std::sqrt(l22Squared);
    constexpr int contourPoints = 121;
    constexpr double deltaChi2 = 2.30;
    const double radius = std::sqrt(deltaChi2);
    contour = std::make_unique<TGraph>(contourPoints);
    contour->SetName(Form("rootfit_contour_%s", dataSet.label));
    for (int i = 0; i < contourPoints; ++i) {
        const double angle = 2. * TMath::Pi() * i / (contourPoints - 1.);
        const double u = radius * std::cos(angle);
        const double v = radius * std::sin(angle);
        contour->SetPoint(i, output.q + l11 * u,
                          output.temp + l21 * u + l22 * v);
    }
    contour->SetTitle(dataSet.label);
    contour->SetLineColor(dataSet.color);
    contour->SetLineWidth(2);
    return true;
}

void qt_contours_rootfit()
{
    constexpr const char *inputPath =
        "HEPData-ins2061074-v1-Figure_16.root";
    constexpr const char *directoryName = "Figure 16";
    constexpr double xMin = 0.;
    constexpr double xMax = 10.;

    std::vector<RootFitDataSet> dataSets = {
        {"Graph1D_y1", "Hist1D_y1_e1plus", "Hist1D_y1_e1minus",
         "Hist1D_y1_e2plus", "Hist1D_y1_e2minus", "0-20%", kRed + 1, {}},
        {"Graph1D_y2", "Hist1D_y2_e1plus", "Hist1D_y2_e1minus",
         "Hist1D_y2_e2plus", "Hist1D_y2_e2minus", "20-40%", kBlue + 1, {}},
        {"Graph1D_y3", "Hist1D_y3_e1plus", "Hist1D_y3_e1minus",
         "Hist1D_y3_e2plus", "Hist1D_y3_e2minus", "40-60%", kGreen + 2, {}},
        {"Graph1D_y4", "Hist1D_y4_e1plus", "Hist1D_y4_e1minus",
         "Hist1D_y4_e2plus", "Hist1D_y4_e2minus", "60-80%", kMagenta + 1, {}},
    };

    const int oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
    if (!input || input->IsZombie()) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot open " << inputPath << '\n';
        return;
    }
    auto *directory =
        dynamic_cast<TDirectoryFile *>(input->Get(directoryName));
    if (!directory) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot find directory \""
                  << directoryName << "\"\n";
        return;
    }
    for (auto &dataSet : dataSets) {
        if (!rootfit_load_dataset(directory, dataSet, xMin, xMax)) {
            gErrorIgnoreLevel = oldErrorIgnoreLevel;
            return;
        }
    }
    input.reset();
    gErrorIgnoreLevel = oldErrorIgnoreLevel;

    std::vector<RootFitResult> fitResults;
    std::vector<std::unique_ptr<TGraph>> contours;
    double previousNor = 100.;
    for (auto &dataSet : dataSets) {
        RootFitResult result;
        std::unique_ptr<TGraph> contour;
        const double norMax = previousNor * (1. - 1.e-6);
        if (!rootfit_dataset(dataSet, xMin, xMax, norMax, result, contour)) {
            return;
        }
        previousNor = result.nor;
        fitResults.push_back(result);
        contours.push_back(std::move(contour));

        std::cout << std::setprecision(10)
                  << dataSet.label
                  << "\nNor = " << result.nor << " +/- " << result.norError
                  << "\nq = " << result.q << " +/- " << result.qError
                  << "\nT = " << result.temp << " +/- " << result.tempError
                  << "\nROOT chi2 = " << result.chi2
                  << "\nndf = " << result.ndf
                  << "\nfit status = " << result.status
                  << "\ncovariance status = " << result.covarianceStatus
                  << "\n\n";
    }

    gStyle->SetOptStat(0);
    auto *canvas = new TCanvas(
        "c_qt_contours_rootfit", "ROOT-fit q-T contours", 900, 700);
    auto *frame = new TH2D(
        "rootfit_frame", "", 10, 1.0, 1.2, 10, 0.0, 0.30);
    frame->GetXaxis()->SetTitle("q");
    frame->GetYaxis()->SetTitle("T (GeV)");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();

    auto *legend = new TLegend(0.68, 0.68, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    for (auto &contour : contours) {
        contour->Draw("L SAME");
        legend->AddEntry(contour.get(), contour->GetTitle(), "l");
    }

    std::vector<std::unique_ptr<TGraph>> bestFitPoints;
    for (std::size_t i = 0; i < fitResults.size(); ++i) {
        auto point = std::make_unique<TGraph>(1);
        point->SetPoint(0, fitResults[i].q, fitResults[i].temp);
        point->SetMarkerStyle(20);
        point->SetMarkerSize(1.2);
        point->SetMarkerColor(dataSets[i].color);
        point->Draw("P SAME");
        bestFitPoints.push_back(std::move(point));
    }
    legend->AddEntry(bestFitPoints.front().get(),
                     "Best-fit points (matching colors)", "p");
    legend->Draw();

    canvas->SaveAs("qt_contours_rootfit_figure16_0_10.png");
    canvas->SaveAs("qt_contours_rootfit_figure16_0_10.pdf");
}
