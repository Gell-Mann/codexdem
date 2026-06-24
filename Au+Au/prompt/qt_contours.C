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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

struct FitPoint {
	double x;
	double y;
	double stat;
	double sys;
};

struct DataSet {
	const char *graphName;
	const char *statName;
	const char *sysName;
	const char *label;
	int color;
	double norMin;
	double norMax;
	std::vector<FitPoint> points;
};

struct ScanResult {
	double nor = 0.;
	double q = 0.;
	double temp = 0.;
	double chi2 = 0.;
	double chi2st = 0.;
	int ndf = 0;
};

double shape_value(double x, double q, double temp)
{
	const double a = 0.;
	const double b = TMath::Sqrt(TMath::Power(x, 2) + TMath::Power(a, 2));
	const double z1 = TMath::Power(1 + (q - 1.) * b / temp, q / (1. - q));
	const double z2 = 2. * b / TMath::Power(2. * TMath::Pi(), 3);
	return z1 * z2 * TMath::Power(10., 4);
}

double best_nor_for_qt(const DataSet &dataSet, double q, double temp)
{
	double numerator = 0.;
	double denominator = 0.;
	for (const auto &point : dataSet.points) {
		const double sigma2 = point.stat * point.stat + point.sys * point.sys;
		const double shape = shape_value(point.x, q, temp);
		numerator += shape * point.y / sigma2;
		denominator += shape * shape / sigma2;
	}

	if (denominator <= 0.) return 0.;
	const double nor = numerator / denominator;
	return TMath::Min(dataSet.norMax, TMath::Max(dataSet.norMin, nor));
}

double total_chi2(const DataSet &dataSet, double nor, double q, double temp)
{
	double chi2 = 0.;
	for (const auto &point : dataSet.points) {
		const double sigma2 = point.stat * point.stat + point.sys * point.sys;
		const double residual = nor * shape_value(point.x, q, temp) - point.y;
		chi2 += residual * residual / sigma2;
	}
	return chi2;
}

double stat_chi2(const DataSet &dataSet, double nor, double q, double temp)
{
	double chi2 = 0.;
	for (const auto &point : dataSet.points) {
		const double residual = nor * shape_value(point.x, q, temp) - point.y;
		chi2 += residual * residual / (point.stat * point.stat);
	}
	return chi2;
}

bool load_dataset(TDirectoryFile *directory, DataSet &dataSet, double xMin, double xMax)
{
	auto *graph = dynamic_cast<TGraphAsymmErrors *>(directory->Get(dataSet.graphName));
	auto *stat = dynamic_cast<TH1F *>(directory->Get(dataSet.statName));
	auto *sys = dynamic_cast<TH1F *>(directory->Get(dataSet.sysName));
	if (!graph || !stat || !sys) {
		std::cerr << "Error: cannot load objects for " << dataSet.label << '\n';
		return false;
	}

	dataSet.points.clear();
	for (int i = 0; i < graph->GetN(); ++i) {
		const double x = graph->GetPointX(i);
		if (x < xMin || x > xMax) continue;
		dataSet.points.push_back({x, graph->GetPointY(i), stat->GetBinError(i + 1), sys->GetBinError(i + 1)});
	}
	return !dataSet.points.empty();
}

std::unique_ptr<TH2D> build_delta_chi2_hist(const DataSet &dataSet, ScanResult &result)
{
	const int qBins = 260;
	const int tBins = 260;
	const double qMin = 1.001;
	const double qMax = 1.5;
	const double tMin = 0.001;
	const double tMax = 0.18;

	auto hist = std::make_unique<TH2D>(Form("delta_%s", dataSet.label), dataSet.label,
					   qBins, qMin, qMax, tBins, tMin, tMax);

	result.chi2 = 1.e300;
	for (int ix = 1; ix <= qBins; ++ix) {
		const double q = hist->GetXaxis()->GetBinCenter(ix);
		for (int iy = 1; iy <= tBins; ++iy) {
			const double temp = hist->GetYaxis()->GetBinCenter(iy);
			const double nor = best_nor_for_qt(dataSet, q, temp);
			const double chi2 = total_chi2(dataSet, nor, q, temp);
			hist->SetBinContent(ix, iy, chi2);
			if (chi2 < result.chi2) {
				result.chi2 = chi2;
				result.nor = nor;
				result.q = q;
				result.temp = temp;
			}
		}
	}

	result.chi2st = stat_chi2(dataSet, result.nor, result.q, result.temp);
	result.ndf = static_cast<int>(dataSet.points.size()) - 3;

	for (int ix = 1; ix <= qBins; ++ix) {
		for (int iy = 1; iy <= tBins; ++iy) {
			hist->SetBinContent(ix, iy, hist->GetBinContent(ix, iy) - result.chi2);
		}
	}

	hist->SetLineColor(dataSet.color);
	hist->SetLineWidth(2);
	hist->SetContour(1);
	double contourLevel[1] = {2.30};
	hist->SetContour(1, contourLevel);
	return hist;
}

void qt_contours()
{
	constexpr const char *inputPath = "HEPData-ins2061074-v1-Figure_15.root";
	constexpr const char *directoryName = "Figure 15";
	constexpr double xMin = 0.;
	constexpr double xMax = 10.;

	std::vector<DataSet> dataSets = {
		{"Graph1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2", "y1", kRed + 1, 0., 100., {}},
		{"Graph1D_y2", "Hist1D_y2_e1", "Hist1D_y2_e2", "y2", kBlue + 1, 3.36271, 12., {}},
		{"Graph1D_y3", "Hist1D_y3_e1", "Hist1D_y3_e2", "y3", kGreen + 2, 1.34511, 9., {}},
		{"Graph1D_y4", "Hist1D_y4_e1", "Hist1D_y4_e2", "y4", kMagenta + 1, 0., 2.5, {}},
	};

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

	std::vector<std::unique_ptr<TH2D>> contourHists;
	std::vector<ScanResult> scanResults;
	for (auto &dataSet : dataSets) {
		ScanResult result;
		contourHists.push_back(build_delta_chi2_hist(dataSet, result));
		scanResults.push_back(result);

		std::cout << dataSet.label << '\n';
		std::cout << "N = " << dataSet.points.size() << '\n';
		std::cout << "Nor = " << result.nor << '\n';
		std::cout << "q = " << result.q << '\n';
		std::cout << "T = " << result.temp << '\n';
		std::cout << "chi2 = " << result.chi2 << '\n';
		std::cout << "chi2st = " << result.chi2st << '\n';
		std::cout << "ndf = " << result.ndf << "\n\n";
	}

	gStyle->SetOptStat(0);
	auto *canvas = new TCanvas("c_qt_contours", "q-T contours", 900, 700);
	auto *frame = new TH2D("frame", "", 10, 1.0, 1.2, 10, 0.0, 0.2);
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

	auto *bestFitPoints = new TGraph(static_cast<int>(scanResults.size()));
	for (std::size_t i = 0; i < scanResults.size(); ++i) {
		bestFitPoints->SetPoint(static_cast<int>(i), scanResults[i].q, scanResults[i].temp);
	}
	bestFitPoints->SetMarkerStyle(20);
	bestFitPoints->SetMarkerSize(1.2);
	bestFitPoints->SetMarkerColor(kRed);
	bestFitPoints->Draw("P SAME");
	legend->AddEntry(bestFitPoints, "Best-fit points", "p");
	legend->Draw();

	canvas->SaveAs("qt_contours_0_10.png");
	canvas->SaveAs("qt_contours_0_10.pdf");
}
