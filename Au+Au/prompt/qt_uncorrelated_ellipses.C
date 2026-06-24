#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct EllipseParameters {
	std::string label;
	double q = 0.;
	double qError = 0.;
	double temp = 0.;
	double tempError = 0.;
	int color = kBlack;
};

bool read_parameters(const char *path, std::vector<EllipseParameters> &parameters)
{
	std::ifstream input(path);
	if (!input) {
		std::cerr << "Error: cannot open " << path << '\n';
		return false;
	}

	const int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
	EllipseParameters *current = nullptr;
	std::string line;
	while (std::getline(input, line)) {
		if (line == "tsa1" || line == "tsa2" || line == "tsa3" || line == "tsa4") {
			parameters.push_back({});
			current = &parameters.back();
			current->label = line;
			current->color = colors[parameters.size() - 1];
			continue;
		}
		if (!current) continue;

		if (std::sscanf(line.c_str(), "q = %lf +/- %lf", &current->q, &current->qError) == 2)
			continue;
		std::sscanf(line.c_str(), "T = %lf +/- %lf", &current->temp, &current->tempError);
	}

	if (parameters.size() != 4) {
		std::cerr << "Error: expected four tsa parameter sets, found " << parameters.size() << '\n';
		return false;
	}
	for (const auto &parameter : parameters) {
		if (parameter.qError <= 0. || parameter.tempError <= 0.) {
			std::cerr << "Error: invalid q or T error for " << parameter.label << '\n';
			return false;
		}
	}
	return true;
}

void qt_uncorrelated_ellipses()
{
	constexpr const char *parameterPath = "tsa_parameters_one_parameter_errors.txt";
	constexpr double deltaChi2 = 2.30;
	constexpr int numberOfPoints = 360;

	std::vector<EllipseParameters> parameters;
	if (!read_parameters(parameterPath, parameters)) return;

	gStyle->SetOptStat(0);
	auto *canvas = new TCanvas("c_qt_uncorrelated", "Uncorrelated q-T ellipses", 900, 700);
	canvas->SetTickx();
	canvas->SetTicky();

	auto *frame = new TH2D("qt_uncorrelated_frame", "", 10, 1.0, 1.2, 10, 0.0, 0.2);
	frame->GetXaxis()->SetTitle("q");
	frame->GetYaxis()->SetTitle("T (GeV)");
	frame->GetXaxis()->CenterTitle();
	frame->GetYaxis()->CenterTitle();
	frame->Draw();

	auto *legend = new TLegend(0.66, 0.66, 0.88, 0.88);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);

	std::vector<TGraph *> ellipses;
	auto *centers = new TGraph(static_cast<int>(parameters.size()));
	const double scale = TMath::Sqrt(deltaChi2);
	for (std::size_t i = 0; i < parameters.size(); ++i) {
		const auto &parameter = parameters[i];
		auto *ellipse = new TGraph(numberOfPoints + 1);
		for (int point = 0; point <= numberOfPoints; ++point) {
			const double angle = 2. * TMath::Pi() * point / numberOfPoints;
			const double q = parameter.q + scale * parameter.qError * TMath::Cos(angle);
			const double temp = parameter.temp + scale * parameter.tempError * TMath::Sin(angle);
			ellipse->SetPoint(point, q, temp);
		}
		ellipse->SetLineColor(parameter.color);
		ellipse->SetLineWidth(2);
		ellipse->Draw("L SAME");
		legend->AddEntry(ellipse, parameter.label.c_str(), "l");
		ellipses.push_back(ellipse);
		centers->SetPoint(static_cast<int>(i), parameter.q, parameter.temp);
	}

	centers->SetMarkerStyle(20);
	centers->SetMarkerSize(1.2);
	centers->SetMarkerColor(kRed);
	centers->Draw("P SAME");
	legend->AddEntry(centers, "Fit points", "p");
	legend->Draw();

	canvas->SaveAs("qt_uncorrelated_ellipses.png");
	canvas->SaveAs("qt_uncorrelated_ellipses.pdf");
}
