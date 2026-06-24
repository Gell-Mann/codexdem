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

struct FitPoint {
    // 单个 pT 数据点及其统计、系统误差。
    double x;
    double y;
    double stat;
    double sys;
};

struct DataSet {
    // Figure 16 中每组数据对应的 ROOT 对象名称和绘图属性。
    const char *graphName;
    const char *statPlusName;
    const char *statMinusName;
    const char *sysPlusName;
    const char *sysMinusName;
    const char *label;
    int color;
    std::vector<FitPoint> points;
};

struct ScanResult {
    // 二维扫描得到的最佳拟合参数及卡方信息。
    double nor = 0.;
    double q = 0.;
    double temp = 0.;
    double chi2 = 0.;
    double chi2st = 0.;
    int ndf = 0;
    bool normalizationConstraintActive = false;
};

double shape_value(double x, double q, double temp)
{
    // 与参考脚本一致的 Tsallis 形状；总归一化由 nor 单独给出。
    const double b = x;
    const double z1 = TMath::Power(1. + (q - 1.) * b / temp, q / (1. - q));
    const double z2 = 2. * b / TMath::Power(2. * TMath::Pi(), 3);
    return z1 * z2 * TMath::Power(10., 4);
}

double fit_variance(const FitPoint &point)
{
    // 统计误差和系统误差按不相关误差进行平方相加。
    return point.stat * point.stat + point.sys * point.sys;
}

double best_nor_for_qt(const DataSet &dataSet, double q, double temp,
                       double norMax = std::numeric_limits<double>::infinity())
{
    // 固定 (q,T) 后，利用加权最小二乘解析求最佳归一化。
    double numerator = 0.;
    double denominator = 0.;
    for (const auto &point : dataSet.points) {
        const double sigma2 = fit_variance(point);
        const double shape = shape_value(point.x, q, temp);
        numerator += shape * point.y / sigma2;
        denominator += shape * shape / sigma2;
    }
    // 物理谱的归一化不允许为负。
    if (denominator <= 0.) return 0.;
    return TMath::Min(norMax, TMath::Max(0., numerator / denominator));
}

double total_chi2(const DataSet &dataSet, double nor, double q, double temp)
{
    // 卡方分母为统计误差和系统误差的平方和。
    double chi2 = 0.;
    for (const auto &point : dataSet.points) {
        const double sigma2 = fit_variance(point);
        const double residual = nor * shape_value(point.x, q, temp) - point.y;
        chi2 += residual * residual / sigma2;
    }
    return chi2;
}

double stat_chi2(const DataSet &dataSet, double nor, double q, double temp)
{
    // 单独计算仅含统计误差的卡方，供结果诊断使用。
    double chi2 = 0.;
    for (const auto &point : dataSet.points) {
        const double residual = nor * shape_value(point.x, q, temp) - point.y;
        chi2 += residual * residual / (point.stat * point.stat);
    }
    return chi2;
}

double symmetric_error(const TH1F *plus, const TH1F *minus, int bin)
{
    // HEPData 分别保存上下误差；这里取绝对值平均得到对称误差。
    return 0.5 * (std::abs(plus->GetBinContent(bin)) +
                  std::abs(minus->GetBinContent(bin)));
}

bool load_dataset(TDirectoryFile *directory, DataSet &dataSet,
                  double xMin, double xMax)
{
    // 读取中心值以及统计/系统误差的上下分量。
    auto *graph = dynamic_cast<TGraphAsymmErrors *>(directory->Get(dataSet.graphName));
    auto *statPlus = dynamic_cast<TH1F *>(directory->Get(dataSet.statPlusName));
    auto *statMinus = dynamic_cast<TH1F *>(directory->Get(dataSet.statMinusName));
    auto *sysPlus = dynamic_cast<TH1F *>(directory->Get(dataSet.sysPlusName));
    auto *sysMinus = dynamic_cast<TH1F *>(directory->Get(dataSet.sysMinusName));
    if (!graph || !statPlus || !statMinus || !sysPlus || !sysMinus) {
        std::cerr << "Error: cannot load objects for " << dataSet.label << '\n';
        return false;
    }

    dataSet.points.clear();
    for (int i = 0; i < graph->GetN(); ++i) {
        const double x = graph->GetPointX(i);
        if (x < xMin || x > xMax) continue;
        const int bin = i + 1;
        const double stat = symmetric_error(statPlus, statMinus, bin);
        const double sys = symmetric_error(sysPlus, sysMinus, bin);
        // 零误差点无法定义卡方权重，因此明确跳过并给出警告。
        if (stat <= 0. || stat * stat + sys * sys <= 0.) {
            std::cerr << "Warning: skip point with zero uncertainty: "
                      << dataSet.label << ", x=" << x << '\n';
            continue;
        }
        dataSet.points.push_back({x, graph->GetPointY(i), stat, sys});
    }
    return !dataSet.points.empty();
}

std::unique_ptr<TH2D> build_delta_chi2_hist(const DataSet &dataSet,
                                             ScanResult &result,
                                             double norMax = std::numeric_limits<double>::infinity())
{
    // q-T 网格范围覆盖 Figure 16 的最佳拟合点及完整置信轮廓。
    constexpr int qBins = 260;
    constexpr int tBins = 430;
    constexpr double qMin = 1.001;
    constexpr double qMax = 1.5;
    constexpr double tMin = 0.001;
    constexpr double tMax = 0.30;

    auto hist = std::make_unique<TH2D>(Form("delta_%s", dataSet.label),
                                       dataSet.label, qBins, qMin, qMax,
                                       tBins, tMin, tMax);
    // 对每个网格点剖面化归一化，只保留 q 和 T 两个轮廓参数。
    result.chi2 = 1.e300;
    for (int ix = 1; ix <= qBins; ++ix) {
        const double q = hist->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= tBins; ++iy) {
            const double temp = hist->GetYaxis()->GetBinCenter(iy);
            const double nor = best_nor_for_qt(dataSet, q, temp, norMax);
            const double chi2 = total_chi2(dataSet, nor, q, temp);
            hist->SetBinContent(ix, iy, chi2);
            if (chi2 < result.chi2) {
                result = {nor, q, temp, chi2, 0., 0};
            }
        }
    }

    result.chi2st = stat_chi2(dataSet, result.nor, result.q, result.temp);
    result.ndf = static_cast<int>(dataSet.points.size()) - 3;
    result.normalizationConstraintActive =
        std::isfinite(norMax) &&
        std::abs(result.nor - norMax) <= 1.e-8 * TMath::Max(1., norMax);
    // 将绝对卡方转换为相对全局最小值的 Delta chi2。
    for (int ix = 1; ix <= qBins; ++ix) {
        for (int iy = 1; iy <= tBins; ++iy) {
            hist->SetBinContent(ix, iy,
                hist->GetBinContent(ix, iy) - result.chi2);
        }
    }
    hist->SetLineColor(dataSet.color);
    hist->SetLineWidth(2);
    // 两个拟合参数的 68.3% 联合置信区域：Delta chi2 = 2.30。
    const double contourLevel[1] = {2.30};
    hist->SetContour(1, contourLevel);
    return hist;
}

void qt_contours()
{
    // 输入文件、ROOT 目录和参与拟合的 pT 范围。
    constexpr const char *inputPath = "HEPData-ins2061074-v1-Figure_16.root";
    constexpr const char *directoryName = "Figure 16";
    constexpr double xMin = 0.;
    constexpr double xMax = 10.;

    // y1--y4 对应从中心到外周的中心度区间。
    std::vector<DataSet> dataSets = {
        {"Graph1D_y1", "Hist1D_y1_e1plus", "Hist1D_y1_e1minus",
         "Hist1D_y1_e2plus", "Hist1D_y1_e2minus", "0-20%", kRed + 1, {}},
        {"Graph1D_y2", "Hist1D_y2_e1plus", "Hist1D_y2_e1minus",
         "Hist1D_y2_e2plus", "Hist1D_y2_e2minus", "20-40%", kBlue + 1, {}},
        {"Graph1D_y3", "Hist1D_y3_e1plus", "Hist1D_y3_e1minus",
         "Hist1D_y3_e2plus", "Hist1D_y3_e2minus", "40-60%", kGreen + 2, {}},
        {"Graph1D_y4", "Hist1D_y4_e1plus", "Hist1D_y4_e1minus",
         "Hist1D_y4_e2plus", "Hist1D_y4_e2minus", "60-80%", kMagenta + 1, {}},
    };

    // HEPData ROOT 文件退出时可能产生对象清理警告，读取期间仅显示致命错误。
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
    for (auto &dataSet : dataSets) {
        if (!load_dataset(directory, dataSet, xMin, xMax)) {
            gErrorIgnoreLevel = oldErrorIgnoreLevel;
            return;
        }
    }
    input.reset();
    gErrorIgnoreLevel = oldErrorIgnoreLevel;

    // 按中心度顺序扫描；后一组 Nor 必须严格小于前一组最佳 Nor。
    std::vector<std::unique_ptr<TH2D>> contourHists;
    std::vector<ScanResult> scanResults;
    double previousNor = std::numeric_limits<double>::infinity();
    for (auto &dataSet : dataSets) {
        ScanResult result;
        const double norMax = std::isfinite(previousNor)
            ? previousNor * (1. - 1.e-6)
            : previousNor;
        contourHists.push_back(build_delta_chi2_hist(dataSet, result, norMax));
        scanResults.push_back(result);
        previousNor = result.nor;
        std::cout << std::setprecision(10)
                  << dataSet.label << "\nN = " << dataSet.points.size()
                  << "\nNor = " << result.nor << "\nq = " << result.q
                  << "\nT = " << result.temp << "\nchi2 = " << result.chi2
                  << "\nchi2st = " << result.chi2st
                  << "\nndf = " << result.ndf
                  << "\nnormalization constraint = "
                  << (result.normalizationConstraintActive ? "active" : "inactive")
                  << "\n\n";
    }

    // 绘制 Delta chi2 轮廓。
    gStyle->SetOptStat(0);
    auto *canvas = new TCanvas("c_qt_contours", "q-T contours", 900, 700);
    auto *frame = new TH2D("frame", "", 10, 1.0, 1.2, 10, 0.0, 0.30);
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

    // 每个最佳拟合点使用其对应轮廓的颜色。
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

    canvas->SaveAs("qt_contours_figure16_0_10.png");
    canvas->SaveAs("qt_contours_figure16_0_10.pdf");
}
