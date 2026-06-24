#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TError.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
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

struct StatFitPoint {
    double x;
    double y;
    double stat;
};

struct StatDataSet {
    const char *graphName;
    const char *statPlusName;
    const char *statMinusName;
    const char *label;
    int color;
    std::vector<StatFitPoint> points;
};

struct StatScanResult {
    double nor = 0.;
    double q = 0.;
    double temp = 0.;
    double chi2 = 0.;
    int ndf = 0;
    bool normalizationConstraintActive = false;
};

double stat_shape_value(double x, double q, double temp)
{
    // 与 qt_contours.C 相同的 Tsallis 形状。
    const double b = x;
    const double z1 = TMath::Power(1. + (q - 1.) * b / temp,
                                   q / (1. - q));
    const double z2 = 2. * b / TMath::Power(2. * TMath::Pi(), 3);
    return z1 * z2 * TMath::Power(10., 4);
}

double stat_best_nor_for_qt(
    const StatDataSet &dataSet, double q, double temp,
    double norMax = std::numeric_limits<double>::infinity())
{
    // 固定 (q,T) 后，以统计误差为权重解析求最佳归一化。
    double numerator = 0.;
    double denominator = 0.;
    for (const auto &point : dataSet.points) {
        const double sigma2 = point.stat * point.stat;
        const double shape = stat_shape_value(point.x, q, temp);
        numerator += shape * point.y / sigma2;
        denominator += shape * shape / sigma2;
    }

    if (denominator <= 0.) return 0.;
    return TMath::Min(norMax, TMath::Max(0., numerator / denominator));
}

double stat_total_chi2(const StatDataSet &dataSet, double nor,
                       double q, double temp)
{
    // 与 qt_contours.C 的核心区别：分母仅为统计误差的平方。
    double chi2 = 0.;
    for (const auto &point : dataSet.points) {
        const double residual =
            nor * stat_shape_value(point.x, q, temp) - point.y;
        chi2 += residual * residual / (point.stat * point.stat);
    }
    return chi2;
}

double stat_symmetric_error(const TH1F *plus, const TH1F *minus, int bin)
{
    // HEPData 分别存放上下误差，取绝对值平均作为对称统计误差。
    return 0.5 * (std::abs(plus->GetBinContent(bin)) +
                  std::abs(minus->GetBinContent(bin)));
}

bool stat_load_dataset(TDirectoryFile *directory, StatDataSet &dataSet,
                       double xMin, double xMax)
{
    auto *graph =
        dynamic_cast<TGraphAsymmErrors *>(directory->Get(dataSet.graphName));
    auto *statPlus =
        dynamic_cast<TH1F *>(directory->Get(dataSet.statPlusName));
    auto *statMinus =
        dynamic_cast<TH1F *>(directory->Get(dataSet.statMinusName));
    if (!graph || !statPlus || !statMinus) {
        std::cerr << "Error: cannot load objects for " << dataSet.label << '\n';
        return false;
    }

    dataSet.points.clear();
    for (int i = 0; i < graph->GetN(); ++i) {
        const double x = graph->GetPointX(i);
        if (x < xMin || x > xMax) continue;

        const double stat = stat_symmetric_error(statPlus, statMinus, i + 1);
        if (stat <= 0.) {
            std::cerr << "Warning: skip point with zero statistical uncertainty: "
                      << dataSet.label << ", x=" << x << '\n';
            continue;
        }
        dataSet.points.push_back({x, graph->GetPointY(i), stat});
    }
    return !dataSet.points.empty();
}

std::unique_ptr<TH2D> stat_build_delta_chi2_hist(
    const StatDataSet &dataSet, StatScanResult &result,
    double norMax = std::numeric_limits<double>::infinity())
{
    // 扫描范围和网格密度与 qt_contours.C 保持一致。
    constexpr int qBins = 260;
    constexpr int tBins = 430;
    constexpr double qMin = 1.001;
    constexpr double qMax = 1.5;
    constexpr double tMin = 0.001;
    constexpr double tMax = 0.30;

    auto hist = std::make_unique<TH2D>(
        Form("stat_delta_%s", dataSet.label), dataSet.label,
        qBins, qMin, qMax, tBins, tMin, tMax);

    result.chi2 = 1.e300;
    for (int ix = 1; ix <= qBins; ++ix) {
        const double q = hist->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= tBins; ++iy) {
            const double temp = hist->GetYaxis()->GetBinCenter(iy);
            const double nor =
                stat_best_nor_for_qt(dataSet, q, temp, norMax);
            const double chi2 = stat_total_chi2(dataSet, nor, q, temp);
            hist->SetBinContent(ix, iy, chi2);
            if (chi2 < result.chi2) {
                result.nor = nor;
                result.q = q;
                result.temp = temp;
                result.chi2 = chi2;
            }
        }
    }

    result.ndf = static_cast<int>(dataSet.points.size()) - 3;
    result.normalizationConstraintActive =
        std::isfinite(norMax) &&
        std::abs(result.nor - norMax) <= 1.e-8 * TMath::Max(1., norMax);

    // 转换为 Delta chi2，并画两个参数的 68.3% 联合置信轮廓。
    for (int ix = 1; ix <= qBins; ++ix) {
        for (int iy = 1; iy <= tBins; ++iy) {
            hist->SetBinContent(
                ix, iy, hist->GetBinContent(ix, iy) - result.chi2);
        }
    }
    hist->SetLineColor(dataSet.color);
    hist->SetLineWidth(2);
    const double contourLevel[1] = {2.30};
    hist->SetContour(1, contourLevel);
    return hist;
}

void qt_contours_stat()
{
    constexpr const char *inputPath =
        "HEPData-ins2061074-v1-Figure_16.root";
    constexpr const char *directoryName = "Figure 16";
    constexpr double xMin = 0.;
    constexpr double xMax = 10.;

    std::vector<StatDataSet> dataSets = {
        {"Graph1D_y1", "Hist1D_y1_e1plus", "Hist1D_y1_e1minus",
         "0-20%", kRed + 1, {}},
        {"Graph1D_y2", "Hist1D_y2_e1plus", "Hist1D_y2_e1minus",
         "20-40%", kBlue + 1, {}},
        {"Graph1D_y3", "Hist1D_y3_e1plus", "Hist1D_y3_e1minus",
         "40-60%", kGreen + 2, {}},
        {"Graph1D_y4", "Hist1D_y4_e1plus", "Hist1D_y4_e1minus",
         "60-80%", kMagenta + 1, {}},
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
        if (!stat_load_dataset(directory, dataSet, xMin, xMax)) {
            gErrorIgnoreLevel = oldErrorIgnoreLevel;
            return;
        }
    }
    input.reset();
    gErrorIgnoreLevel = oldErrorIgnoreLevel;

    std::vector<std::unique_ptr<TH2D>> contourHists;
    std::vector<StatScanResult> scanResults;
    double previousNor = std::numeric_limits<double>::infinity();
    for (auto &dataSet : dataSets) {
        StatScanResult result;
        const double norMax = std::isfinite(previousNor)
            ? previousNor * (1. - 1.e-6)
            : previousNor;
        contourHists.push_back(
            stat_build_delta_chi2_hist(dataSet, result, norMax));
        scanResults.push_back(result);
        previousNor = result.nor;

        std::cout << std::setprecision(10)
                  << dataSet.label << "\nN = " << dataSet.points.size()
                  << "\nNor = " << result.nor
                  << "\nq = " << result.q
                  << "\nT = " << result.temp
                  << "\nchi2(stat only) = " << result.chi2
                  << "\nndf = " << result.ndf
                  << "\nnormalization constraint = "
                  << (result.normalizationConstraintActive
                          ? "active" : "inactive")
                  << "\n\n";
    }

    gStyle->SetOptStat(0);
    auto *canvas =
        new TCanvas("c_qt_contours_stat", "q-T stat-only contours", 900, 700);
    auto *frame =
        new TH2D("stat_frame", "", 10, 1.0, 1.2, 10, 0.0, 0.30);
    frame->GetXaxis()->SetTitle("q");
    frame->GetYaxis()->SetTitle("T (GeV)");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();

    auto *legend = new TLegend(0.68, 0.68, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    for (auto &hist : contourHists) {
        hist->Draw("CONT3 SAME");
        legend->AddEntry(hist.get(), hist->GetTitle(), "l");
    }

    std::vector<std::unique_ptr<TGraph>> bestFitPoints;
    for (std::size_t i = 0; i < scanResults.size(); ++i) {
        auto point = std::make_unique<TGraph>(1);
        point->SetPoint(0, scanResults[i].q, scanResults[i].temp);
        point->SetMarkerStyle(20);
        point->SetMarkerSize(1.2);
        point->SetMarkerColor(dataSets[i].color);
        point->Draw("P SAME");
        bestFitPoints.push_back(std::move(point));
    }
    legend->AddEntry(bestFitPoints.front().get(),
                     "Best-fit points (matching colors)", "p");
    legend->Draw();

    canvas->SaveAs("qt_contours_stat_figure16_0_10.png");
    canvas->SaveAs("qt_contours_stat_figure16_0_10.pdf");
}
