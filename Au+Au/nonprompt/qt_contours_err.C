#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TError.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
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
    // 单参数 1 sigma 剖面误差（Delta chi2 = 1），上下不对称。
    // 新成员追加在末尾，以保持扫描循环里的位置初始化 {nor,q,temp,chi2,0.,0} 不受影响。
    double qErrLow = 0.;
    double qErrHigh = 0.;
    double tempErrLow = 0.;
    double tempErrHigh = 0.;
    double norErrLow = 0.;
    double norErrHigh = 0.;
    bool qErrValid = false;
    bool tempErrValid = false;
    bool norErrValid = false;
};

// 扫描网格常数提到文件作用域，确保直方图构建与 Nor 误差扫描使用完全相同的网格。
namespace grid {
constexpr int qBins = 260;
constexpr int tBins = 430;
constexpr double qMin = 1.001;
constexpr double qMax = 1.5;
constexpr double tMin = 0.001;
constexpr double tMax = 0.30;
}  // namespace grid

// 在一维剖面曲线 (xs, dchi2) 上寻找 dchi2 = level 的左右交点。
// dchi2 应已减去其最小值；errLow/errHigh 返回相对最小值点的正向距离。
// 若某一侧在扫描区间内未达到 level，则对应 valid 置为 false（说明区间被截断）。
void find_crossings(const std::vector<double> &xs,
                    const std::vector<double> &dchi2, double level,
                    double &errLow, double &errHigh, bool &lowValid,
                    bool &highValid)
{
    errLow = 0.;
    errHigh = 0.;
    lowValid = false;
    highValid = false;
    if (xs.size() < 2) return;

    std::size_t imin = 0;
    for (std::size_t i = 1; i < dchi2.size(); ++i) {
        if (dchi2[i] < dchi2[imin]) imin = i;
    }
    const double best = xs[imin];

    // 向左扫描：找到第一个 dchi2 >= level 的格点，并与其右邻线性插值。
    for (std::size_t k = imin; k > 0;) {
        --k;
        if (dchi2[k] >= level) {
            const double x = xs[k] + (level - dchi2[k]) *
                                         (xs[k + 1] - xs[k]) /
                                         (dchi2[k + 1] - dchi2[k]);
            errLow = best - x;
            lowValid = true;
            break;
        }
    }
    // 向右扫描：找到第一个 dchi2 >= level 的格点，并与其左邻线性插值。
    for (std::size_t k = imin + 1; k < dchi2.size(); ++k) {
        if (dchi2[k] >= level) {
            const double x = xs[k - 1] + (level - dchi2[k - 1]) *
                                             (xs[k] - xs[k - 1]) /
                                             (dchi2[k] - dchi2[k - 1]);
            errHigh = x - best;
            highValid = true;
            break;
        }
    }
}

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
    using grid::qBins;
    using grid::qMax;
    using grid::qMin;
    using grid::tBins;
    using grid::tMax;
    using grid::tMin;

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

    // 单参数误差：在 Nor 已被剖面化的 Delta chi2 曲面上再把另一参数也剖面掉
    // （逐列/逐行取最小），得到一维剖面，找 Delta chi2 = 1 的左右交点。
    constexpr double singleParamLevel = 1.0;
    std::vector<double> qVals(qBins), qProfile(qBins);
    for (int ix = 1; ix <= qBins; ++ix) {
        double best = 1.e300;
        for (int iy = 1; iy <= tBins; ++iy) {
            best = TMath::Min(best, hist->GetBinContent(ix, iy));
        }
        qVals[ix - 1] = hist->GetXaxis()->GetBinCenter(ix);
        qProfile[ix - 1] = best;
    }
    {
        bool lowValid = false, highValid = false;
        find_crossings(qVals, qProfile, singleParamLevel, result.qErrLow,
                       result.qErrHigh, lowValid, highValid);
        result.qErrValid = lowValid && highValid;
    }

    std::vector<double> tVals(tBins), tProfile(tBins);
    for (int iy = 1; iy <= tBins; ++iy) {
        double best = 1.e300;
        for (int ix = 1; ix <= qBins; ++ix) {
            best = TMath::Min(best, hist->GetBinContent(ix, iy));
        }
        tVals[iy - 1] = hist->GetYaxis()->GetBinCenter(iy);
        tProfile[iy - 1] = best;
    }
    {
        bool lowValid = false, highValid = false;
        find_crossings(tVals, tProfile, singleParamLevel, result.tempErrLow,
                       result.tempErrHigh, lowValid, highValid);
        result.tempErrValid = lowValid && highValid;
    }

    hist->SetLineColor(dataSet.color);
    hist->SetLineWidth(2);
    // 两个拟合参数的 68.3% 联合置信区域：Delta chi2 = 2.30。
    const double contourLevel[1] = {2.30};
    hist->SetContour(1, contourLevel);
    return hist;
}

void compute_nor_error(const DataSet &dataSet, ScanResult &result)
{
    // 方案一：条件误差。把 q、T 固定在最佳拟合值 q*、T*，只让 Nor 变化。
    // 此时卡方对 Nor 是严格的抛物线
    //   chi2(Nor) = A*Nor^2 - 2*B*Nor + C,  A = sum S^2/sigma2,
    // 其 Delta chi2 = 1 的半宽解析可得：sigma_Nor = 1 / sqrt(A)。
    //
    // 与剖面误差不同，这里不让 q、T 浮动，因此不受 Nor-q-T 简并影响，
    // 始终给出有限、对称、良定义的归一化不确定度（即把 q、T 视为已知时
    // Nor 的统计误差，等价于 _rootfit 版 Minuit 的对称 ParError）。
    double a = 0.;
    for (const auto &point : dataSet.points) {
        const double sigma2 = fit_variance(point);
        const double shape = shape_value(point.x, result.q, result.temp);
        a += shape * shape / sigma2;
    }
    if (a > 0.) {
        const double sigmaNor = 1. / std::sqrt(a);
        result.norErrLow = sigmaNor;
        result.norErrHigh = sigmaNor;
        result.norErrValid = true;
    } else {
        result.norErrLow = 0.;
        result.norErrHigh = 0.;
        result.norErrValid = false;
    }
}

void qt_contours_err()
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
        // q、T 误差已由曲面剖面给出；Nor 误差需额外的外层扫描。
        compute_nor_error(dataSet, result);
        scanResults.push_back(result);
        previousNor = result.nor;

        // 把 +x/-x 不对称误差格式化输出；若某侧被扫描边界截断则标注。
        auto fmtErr = [](double low, double high, bool valid) {
            std::ostringstream os;
            os << std::setprecision(6) << " (+" << high << " / -" << low << ")";
            if (!valid) os << " [truncated: widen scan range]";
            return os.str();
        };

        std::cout << std::setprecision(10)
                  << dataSet.label << "\nN = " << dataSet.points.size()
                  << "\nNor = " << result.nor
                  << fmtErr(result.norErrLow, result.norErrHigh, result.norErrValid)
                  << "\nq = " << result.q
                  << fmtErr(result.qErrLow, result.qErrHigh, result.qErrValid)
                  << "\nT = " << result.temp
                  << fmtErr(result.tempErrLow, result.tempErrHigh, result.tempErrValid)
                  << "\nchi2 = " << result.chi2
                  << "\nchi2st = " << result.chi2st
                  << "\nndf = " << result.ndf
                  << "\nnormalization constraint = "
                  << (result.normalizationConstraintActive ? "active" : "inactive")
                  << "\n\n";
    }

    const char *resultPath = "qt_contours_err_results.txt";
    std::ofstream resultFile(resultPath);
    if (!resultFile) {
        std::cerr << "Error: cannot write " << resultPath << '\n';
        return;
    }
    resultFile << "# q and T errors: single-parameter 1 sigma "
               << "(Delta chi2 = 1 profile)\n"
               << "# Nor error: 1 sigma at fixed best q,T\n"
               << "# centrality\tNor\tNor_err_low\tNor_err_high"
               << "\tq\tq_err_low\tq_err_high"
               << "\tT\tT_err_low\tT_err_high\n";
    resultFile << std::setprecision(10);
    for (std::size_t i = 0; i < scanResults.size(); ++i) {
        const ScanResult &r = scanResults[i];
        resultFile << dataSets[i].label
                   << '\t' << r.nor
                   << '\t' << r.norErrLow
                   << '\t' << r.norErrHigh
                   << '\t' << r.q
                   << '\t' << r.qErrLow
                   << '\t' << r.qErrHigh
                   << '\t' << r.temp
                   << '\t' << r.tempErrLow
                   << '\t' << r.tempErrHigh
                   << '\n';
    }
    std::cout << "Saved fit parameters to " << resultPath << "\n\n";

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

    // 每个最佳拟合点使用其对应轮廓的颜色，并带上 q、T 的不对称 1σ 误差棒。
    std::vector<std::unique_ptr<TGraphAsymmErrors>> bestFitPoints;
    for (std::size_t i = 0; i < scanResults.size(); ++i) {
        const ScanResult &r = scanResults[i];
        auto point = std::make_unique<TGraphAsymmErrors>(1);
        point->SetPoint(0, r.q, r.temp);
        // x 方向 = q 误差，y 方向 = T 误差（下，上）。
        point->SetPointError(0, r.qErrLow, r.qErrHigh, r.tempErrLow,
                             r.tempErrHigh);
        point->SetMarkerStyle(20);
        point->SetMarkerSize(1.2);
        point->SetMarkerColor(dataSets[i].color);
        point->SetLineColor(dataSets[i].color);
        point->SetLineWidth(2);
        point->Draw("P SAME");
        bestFitPoints.push_back(std::move(point));
    }
    legend->AddEntry(bestFitPoints.front().get(),
                     "Best-fit points (#pm1#sigma, matching colors)", "pl");
    legend->Draw();

    // 在左上空白区用 TLatex 旁注每组的 q、T 及其不对称 1σ 误差，颜色与轮廓一致。
    auto fmtAsym = [](double v, double low, double high, int prec) {
        return std::string(Form("%.*f^{+%.*f}_{-%.*f}", prec, v, prec, high,
                                prec, low));
    };
    auto *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.022);
    double yText = 0.88;
    for (std::size_t i = 0; i < scanResults.size(); ++i) {
        const ScanResult &r = scanResults[i];
        const std::string qStr = fmtAsym(r.q, r.qErrLow, r.qErrHigh, 3);
        const std::string tStr = fmtAsym(r.temp, r.tempErrLow, r.tempErrHigh, 3);
        const std::string line =
            Form("#color[%d]{%s:  q=%s,  T=%s}", dataSets[i].color,
                 dataSets[i].label, qStr.c_str(), tStr.c_str());
        latex->DrawLatex(0.135, yText, line.c_str());
        yText -= 0.045;
    }
    // 说明误差定义，避免误读（Nor 的 1σ 见终端输出）。
    latex->SetTextColor(kGray + 2);
    latex->DrawLatex(0.135, yText - 0.01,
                     "q, T: single-param 1#sigma (#Delta#chi^{2}=1 profile)");
    latex->DrawLatex(0.135, yText - 0.05,
                     "Nor 1#sigma at fixed best q,T: see printout");

    auto *panel = new TLatex();
    panel->SetNDC();
    panel->SetTextFont(42);
    panel->SetTextSize(0.04);
    panel->DrawLatex(0.15, 0.15, "(d)");

    canvas->SaveAs("qt_contours_err_figure16_0_10.png");
    canvas->SaveAs("qt_contours_err_figure16_0_10.pdf");
}
