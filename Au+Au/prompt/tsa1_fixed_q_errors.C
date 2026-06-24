#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TError.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TMath.h>

#include <iostream>
#include <memory>
#include <vector>

struct FitConfig {
	const char *label;
	const char *graphName;
	const char *statErrorName;
	const char *sysErrorName;
	double nor0;
	double norMin;
	double norMax;
	double q0;
	double qMin;
	double qMax;
	double t0;
	double tMin;
	double tMax;
	bool fixQInFullFit;
	bool fixTInFullFit;
};

struct FitObjects {
	TGraphAsymmErrors *gr = nullptr;
	TH1F *e1 = nullptr;
	TH1F *e2 = nullptr;
};

Double_t func1_single_parameter_error(Double_t *x, Double_t *par)
{
	Double_t b = TMath::Sqrt(TMath::Power(x[0], 2));
	Double_t z1 = TMath::Power(1 + (par[1] - 1) * b / par[2], par[1] / (1 - par[1]));
	Double_t z2 = 2 * par[0] * b / TMath::Power(2 * TMath::Pi(), 3);
	return z1 * z2 * TMath::Power(10, 4);
}

void apply_full_fit_settings(TF1 *fit, const FitConfig &config)
{
	fit->SetParNames("Nor", "q", "T");
	fit->SetParameter(0, config.nor0);
	fit->SetParLimits(0, config.norMin, config.norMax);
	if (config.fixQInFullFit) {
		fit->FixParameter(1, config.q0);
	} else {
		fit->SetParameter(1, config.q0);
		fit->SetParLimits(1, config.qMin, config.qMax);
	}
	if (config.fixTInFullFit) {
		fit->FixParameter(2, config.t0);
	} else {
		fit->SetParameter(2, config.t0);
		fit->SetParLimits(2, config.tMin, config.tMax);
	}
}

FitObjects load_fit_objects(TDirectoryFile *directory, const FitConfig &config)
{
	FitObjects objects;

	auto *inputGraph = dynamic_cast<TGraphAsymmErrors *>(directory->Get(config.graphName));
	auto *inputE1 = dynamic_cast<TH1F *>(directory->Get(config.statErrorName));
	auto *inputE2 = dynamic_cast<TH1F *>(directory->Get(config.sysErrorName));
	if (!inputGraph || !inputE1 || !inputE2) {
		std::cerr << "Error: cannot find objects for " << config.label << '\n';
		return objects;
	}

	objects.gr = dynamic_cast<TGraphAsymmErrors *>(inputGraph->Clone(Form("%s_graph", config.label)));
	objects.e1 = dynamic_cast<TH1F *>(inputE1->Clone(Form("%s_stat", config.label)));
	objects.e2 = dynamic_cast<TH1F *>(inputE2->Clone(Form("%s_sys", config.label)));
	if (!objects.gr || !objects.e1 || !objects.e2) {
		std::cerr << "Error: cannot clone objects for " << config.label << '\n';
		return FitObjects{};
	}

	objects.e1->SetDirectory(nullptr);
	objects.e2->SetDirectory(nullptr);
	return objects;
}

void calculate_chi2(TGraphAsymmErrors *gr, TH1F *e1, TH1F *e2, TF1 *fit,
		    double xMin, double xMax, int parN,
		    int &n, int &ndf, double &chi2, double &chi2st)
{
	n = 0;
	chi2 = 0.;
	chi2st = 0.;

	for (int i = 1; i < gr->GetN() + 1; i++) {
		double x = gr->GetPointX(i - 1);
		if (x < xMin || x > xMax) continue;

		n++;
		double residual = fit->Eval(x) - gr->GetPointY(i - 1);
		chi2 += residual * residual / (pow(e1->GetBinError(i), 2) + pow(e2->GetBinError(i), 2));
		chi2st += residual * residual / pow(e1->GetBinError(i), 2);
	}

	ndf = n - parN;
}

double one_parameter_error(TGraphAsymmErrors *gr, const char *name, double xMin, double xMax,
			   double nor, double q, double temp, const FitConfig &config, int freeIndex)
{
	auto *fit = new TF1(name, func1_single_parameter_error, xMin, xMax, 3);
	fit->SetParNames("Nor", "q", "T");

	if (freeIndex == 0) {
		fit->SetParameter(0, nor);
		fit->SetParLimits(0, config.norMin, config.norMax);
	} else {
		fit->FixParameter(0, nor);
	}

	if (freeIndex == 1) {
		fit->SetParameter(1, q);
		fit->SetParLimits(1, config.qMin, config.qMax);
	} else {
		fit->FixParameter(1, q);
	}

	if (freeIndex == 2) {
		fit->SetParameter(2, temp);
		fit->SetParLimits(2, config.tMin, config.tMax);
	} else {
		fit->FixParameter(2, temp);
	}

	gr->Fit(fit, "RQ");
	return fit->GetParError(freeIndex);
}

void tsa1_fixed_q_errors()
{
	constexpr const char *inputPath = "HEPData-ins2061074-v1-Figure_15.root";
	constexpr const char *directoryName = "Figure 15";
	constexpr double xMin = 0.;
	constexpr double xMax = 10.;
	constexpr int parN = 3;

	std::vector<FitConfig> configs = {
		{"tsa1", "Graph1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2",
		 8.363, 0., 100., 1.110, 1.001, 1.5, 0.107, 0.001, 0.18, false, false},
		{"tsa2", "Graph1D_y2", "Hist1D_y2_e1", "Hist1D_y2_e2",
		 6.55588, 3.36271, 12., 1.117, 1.001, 1.5, 0.078, 0.001, 0.18, true, false},
		{"tsa3", "Graph1D_y3", "Hist1D_y3_e1", "Hist1D_y3_e2",
		 2.13295, 1.34511, 9., 1.118, 1.001, 1.5, 0.059, 0.001, 0.18, true, false},
		{"tsa4", "Graph1D_y4", "Hist1D_y4_e1", "Hist1D_y4_e2",
		 0.312712, 0., 2.5, 1.142, 1.001, 1.5, 0.073, 0.001, 0.18, false, true},
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

	auto *canvas = new TCanvas("c_single_parameter_errors", "single parameter errors", 900, 700);
	canvas->SetLogy();

	for (auto &config : configs) {
		FitObjects objects = load_fit_objects(directory, config);
		if (!objects.gr || !objects.e1 || !objects.e2) return;

		auto *fullFit = new TF1(Form("%s_fullFit", config.label), func1_single_parameter_error,
					xMin, xMax, parN);
		apply_full_fit_settings(fullFit, config);

		objects.gr->Fit(fullFit, "RQ");

		double nor = fullFit->GetParameter(0);
		double q = fullFit->GetParameter(1);
		double temp = fullFit->GetParameter(2);

		double norFullErr = fullFit->GetParError(0);
		double qFullErr = fullFit->GetParError(1);
		double tempFullErr = fullFit->GetParError(2);

		double norOneErr = one_parameter_error(objects.gr, Form("%s_norErrorFit", config.label),
						       xMin, xMax, nor, q, temp, config, 0);
		double qOneErr = one_parameter_error(objects.gr, Form("%s_qErrorFit", config.label),
						     xMin, xMax, nor, q, temp, config, 1);
		double tempOneErr = one_parameter_error(objects.gr, Form("%s_tErrorFit", config.label),
							xMin, xMax, nor, q, temp, config, 2);

		int n = 0;
		int ndf = 0;
		double chi2 = 0.;
		double chi2st = 0.;
		calculate_chi2(objects.gr, objects.e1, objects.e2, fullFit, xMin, xMax,
			       parN, n, ndf, chi2, chi2st);

		std::cout << "\n" << config.label << '\n';
		std::cout << "Original full fit correlated errors\n";
		std::cout << "Nor = " << nor << " +/- " << norFullErr << '\n';
		std::cout << "q = " << q << " +/- " << qFullErr << '\n';
		std::cout << "T = " << temp << " +/- " << tempFullErr << '\n';
		std::cout << "N = " << n << '\n';
		std::cout << "chi2 = " << chi2 << '\n';
		std::cout << "chi2st = " << chi2st << '\n';
		std::cout << "ndf = " << ndf << '\n';

		std::cout << "Report values with one-parameter errors\n";
		std::cout << "Nor = " << nor << " +/- " << norOneErr << '\n';
		std::cout << "q = " << q << " +/- " << qOneErr << '\n';
		std::cout << "T = " << temp << " +/- " << tempOneErr << '\n';
	}

	input.reset();
	gErrorIgnoreLevel = oldErrorIgnoreLevel;
}
