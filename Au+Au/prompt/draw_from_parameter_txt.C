#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TDirectoryFile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TError.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace TMath;

struct TsaParameter {
	double nor = 0.;
	double norErr = 0.;
	double q = 0.;
	double qErr = 0.;
	double temp = 0.;
	double tempErr = 0.;
};

Double_t func_from_txt(Double_t *x, Double_t *par)
{
	Double_t a = 0;
	Double_t b = TMath::Sqrt(TMath::Power(x[0], 2) + TMath::Power(a, 2));
	Double_t z1 = TMath::Power(1 + (par[1] - 1) * b / par[2], par[1] / (1 - par[1]));
	Double_t z2 = 2 * par[0] * b / TMath::Power(2 * TMath::Pi(), 3);
	return z1 * z2 * TMath::Power(10, 4);
}

bool read_parameters(const char *path, TsaParameter par[4])
{
	std::ifstream input(path);
	if (!input) {
		std::cerr << "Error: cannot open " << path << '\n';
		return false;
	}

	int index = -1;
	std::string line;
	while (std::getline(input, line)) {
		if (line == "tsa1") index = 0;
		else if (line == "tsa2") index = 1;
		else if (line == "tsa3") index = 2;
		else if (line == "tsa4") index = 3;
		else if (index >= 0) {
			std::istringstream stream(line);
			std::string name;
			std::string eq;
			std::string pm;
			double value = 0.;
			double error = 0.;
			stream >> name >> eq >> value >> pm >> error;
			if (!stream || eq != "=" || pm != "+/-") continue;
			if (name == "Nor") {
				par[index].nor = value;
				par[index].norErr = error;
			} else if (name == "q") {
				par[index].q = value;
				par[index].qErr = error;
			} else if (name == "T") {
				par[index].temp = value;
				par[index].tempErr = error;
			}
		}
	}

	for (int i = 0; i < 4; ++i) {
		if (par[i].nor <= 0. || par[i].q <= 0. || par[i].temp <= 0.) {
			std::cerr << "Error: incomplete parameters for tsa" << i + 1 << '\n';
			return false;
		}
	}
	return true;
}

TGraphAsymmErrors *clone_graph(TDirectoryFile *td, const char *name, const char *cloneName)
{
	auto *input = dynamic_cast<TGraphAsymmErrors *>(td->Get(name));
	if (!input) {
		std::cerr << "Error: cannot find " << name << '\n';
		return nullptr;
	}
	return dynamic_cast<TGraphAsymmErrors *>(input->Clone(cloneName));
}

void scale_graph(TGraphAsymmErrors *graph, double scale)
{
	for (int i = 0; i < graph->GetN(); i++) {
		double x, y;
		graph->GetPoint(i, x, y);
		graph->SetPoint(i, x, scale * y);
		graph->SetPointEYlow(i, scale * graph->GetErrorYlow(i));
		graph->SetPointEYhigh(i, scale * graph->GetErrorYhigh(i));
	}
}

TF1 *make_fixed_curve(const char *name, double scale, const TsaParameter &par,
		      double xMin, double xMax)
{
	TF1 *fit = new TF1(name, func_from_txt, xMin, xMax, 3);
	fit->SetLineColor(2);
	fit->SetParNames("Nor", "q", "T");
	fit->FixParameter(0, scale * par.nor);
	fit->FixParameter(1, par.q);
	fit->FixParameter(2, par.temp);
	return fit;
}

void draw_cenm1_from_txt(const TsaParameter par[4], TDirectoryFile *td)
{
	TGraphAsymmErrors *gr1 = clone_graph(td, "Graph1D_y1", "gr1_txt");
	TGraphAsymmErrors *gr2 = clone_graph(td, "Graph1D_y2", "gr2_txt");
	TGraphAsymmErrors *gr3 = clone_graph(td, "Graph1D_y3", "gr3_txt");
	TGraphAsymmErrors *gr4 = clone_graph(td, "Graph1D_y4", "gr4_txt");
	if (!gr1 || !gr2 || !gr3 || !gr4) return;

	TCanvas *c1 = new TCanvas("fit1_txt", "fit1_txt", 960, 0, 550, 500);

	double xMin = 0.;
	double xMax = 12.;

	TF1 *f1 = make_fixed_curve("f1_txt", 30., par[0], xMin, xMax);
	scale_graph(gr1, 30.);

	c1->SetTickx();
	c1->SetTicky();
	c1->SetLogy();

	gr1->GetXaxis()->SetLimits(0, 11);
	gr1->SetMinimum(pow(10, -9));
	gr1->SetMaximum(pow(10, 3));

	gr1->SetLineColor(kCyan - 2);
	gr1->SetLineWidth(2);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(4);
	gr1->SetMarkerSize(1);

	gStyle->SetTitleY(0.96);
	gStyle->SetLineScalePS(2);
	gStyle->SetLineStyleString(10, "15 15");

	gr1->SetTitle("");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr1->GetXaxis()->SetTitleOffset(1.3);
	gr1->GetYaxis()->CenterTitle();
	gr1->GetYaxis()->SetLabelSize(0.03);
	gr1->GetYaxis()->SetTitle("1/N 1/(2#pip_{T}) d^{2}N/(dydp_{T}) ((GeV/c)^{-2})");
	gr1->GetYaxis()->SetTitleOffset(1.3);

	gr1->Draw("AP");
	f1->Draw("same");

	TF1 *f2 = make_fixed_curve("f2_txt", 10., par[1], xMin, xMax);
	scale_graph(gr2, 10.);
	gr2->SetLineColor(kCyan - 2);
	gr2->SetLineWidth(2);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(1);
	gr2->SetMarkerSize(1);
	gr2->Draw("P same");
	f2->Draw("same");

	TF1 *f3 = make_fixed_curve("f3_txt", 5., par[2], xMin, xMax);
	scale_graph(gr3, 5.);
	gr3->SetLineColor(kCyan - 2);
	gr3->SetLineWidth(2);
	gr3->SetMarkerStyle(22);
	gr3->SetMarkerColor(4);
	gr3->SetMarkerSize(1);
	gr3->Draw("P same");
	f3->Draw("same");

	TF1 *f4 = make_fixed_curve("f4_txt", 1., par[3], xMin, xMax);
	scale_graph(gr4, 1.);
	gr4->SetLineColor(kCyan - 2);
	gr4->SetLineWidth(2);
	gr4->SetMarkerStyle(23);
	gr4->SetMarkerColor(1);
	gr4->SetMarkerSize(1);
	gr4->Draw("P same");
	f4->Draw("same");

	TLine *line1 = new TLine(0.20, 0.78, 0.25, 0.78); line1->SetNDC();
	TLatex *tex1 = new TLatex; tex1->SetNDC();
	tex1->SetTextSize(0.03);
	tex1->SetTextFont(42);
	tex1->DrawLatex(0.27, 0.77, "Tsallis");
	line1->SetLineStyle(1);
	line1->SetLineColor(2);
	line1->SetLineWidth(2);
	line1->Draw("same");

	TLatex *tex20 = new TLatex; tex20->SetNDC();
	tex20->SetTextSize(0.03);
	tex20->SetTextFont(132);
	tex20->DrawLatex(0.2, 0.84, "direct-photon 200 GeV Au+Au, PHENIX");

	TLine *line21 = new TLine(0.65, 0.85, 0.70, 0.85); line21->SetNDC();
	line21->SetLineColor(kCyan - 2);
	line21->SetLineWidth(2);
	line21->Draw();
	TMarker *ma = new TMarker((0.65 + 0.70) / 2, 0.85, 20); ma->SetNDC();
	ma->SetMarkerColor(4);
	ma->SetMarkerSize(1);
	ma->Draw("same");
	TLatex *tex21 = new TLatex; tex21->SetNDC();
	tex21->SetTextSize(0.03);
	tex21->SetTextFont(42);
	tex21->DrawLatex(0.71, 0.84, "0#minus20\%, #times30");

	TLine *line22 = new TLine(0.65, 0.80, 0.70, 0.80); line22->SetNDC();
	line22->SetLineColor(kCyan - 2);
	line22->SetLineWidth(2);
	line22->Draw();
	TMarker *ma2 = new TMarker((0.65 + 0.70) / 2, 0.80, 21); ma2->SetNDC();
	ma2->SetMarkerColor(1);
	ma2->SetMarkerSize(1);
	ma2->Draw("same");
	TLatex *tex22 = new TLatex; tex22->SetNDC();
	tex22->SetTextSize(0.03);
	tex22->SetTextFont(42);
	tex22->DrawLatex(0.71, 0.79, "20#minus40\%, #times10");

	TLine *line23 = new TLine(0.65, 0.75, 0.70, 0.75); line23->SetNDC();
	line23->SetLineColor(kCyan - 2);
	line23->SetLineWidth(2);
	line23->Draw();
	TMarker *ma3 = new TMarker((0.65 + 0.70) / 2, 0.75, 22); ma3->SetNDC();
	ma3->SetMarkerColor(4);
	ma3->SetMarkerSize(1);
	ma3->Draw("same");
	TLatex *tex23 = new TLatex; tex23->SetNDC();
	tex23->SetTextSize(0.03);
	tex23->SetTextFont(42);
	tex23->DrawLatex(0.71, 0.74, "40#minus60\%, #times5");

	TLine *line24 = new TLine(0.65, 0.70, 0.70, 0.70); line24->SetNDC();
	line24->SetLineColor(kCyan - 2);
	line24->SetLineWidth(2);
	line24->Draw();
	TMarker *ma4 = new TMarker((0.65 + 0.70) / 2, 0.70, 23); ma4->SetNDC();
	ma4->SetMarkerColor(1);
	ma4->SetMarkerSize(1);
	ma4->Draw("same");
	TLatex *tex24 = new TLatex; tex24->SetNDC();
	tex24->SetTextSize(0.03);
	tex24->SetTextFont(42);
	tex24->DrawLatex(0.71, 0.69, "60#minus93\%, #times1");

	TLatex *tex27 = new TLatex; tex27->SetNDC();
	tex27->SetTextSize(0.04);
	tex27->SetTextFont(42);
	tex27->DrawLatex(0.15, 0.15, "(a)");

	c1->SaveAs("fig1a_from_txt.eps");
}

void draw_q_from_txt(const TsaParameter par[4])
{
	double x1[] = {10, 30, 50, 76.5};
	double y1[] = {par[0].q, par[1].q, par[2].q, par[3].q};
	double ex1[] = {0., 0., 0., 0.};
	double ey1[] = {par[0].qErr, par[1].qErr, par[2].qErr, par[3].qErr};
	int numpoints = 4;

	TGraphErrors *graph1 = new TGraphErrors(numpoints, x1, y1, ex1, ey1);
	TCanvas *c1 = new TCanvas("sca_q_txt", "sca_q_txt", 960, 570, 500, 450);

	gStyle->SetTitleY(0.96);
	graph1->SetTitle("");
	graph1->GetXaxis()->CenterTitle();
	graph1->GetXaxis()->SetTitle("Centrality");
	graph1->GetXaxis()->SetTitleOffset(1.3);
	graph1->GetYaxis()->CenterTitle();
	graph1->GetYaxis()->SetLabelSize(0.03);
	graph1->GetYaxis()->SetTitle("q");
	graph1->GetYaxis()->SetTitleOffset(1.48);
	graph1->SetMarkerColor(4);
	graph1->SetMarkerStyle(20);
	graph1->SetMarkerSize(1.3);
	graph1->Draw("AP");
	c1->SetTickx();
	c1->SetTicky();

	TLatex *tex1 = new TLatex; tex1->SetNDC();
	tex1->SetTextSize(0.04);
	tex1->SetTextFont(42);
	tex1->DrawLatex(0.17, 0.84, "200 GeV Au+Au");
	TLatex *tex2 = new TLatex; tex2->SetNDC();
	tex2->SetTextSize(0.04);
	tex2->SetTextFont(42);
	tex2->DrawLatex(0.83, 0.15, "(b)");

	c1->SaveAs("fig1b_from_txt.eps");
}

void draw_t_from_txt(const TsaParameter par[4])
{
	double x1[] = {10, 30, 50, 76.5};
	double y1[] = {par[0].temp, par[1].temp, par[2].temp, par[3].temp};
	double ex1[] = {0., 0., 0., 0.};
	double ey1[] = {par[0].tempErr, par[1].tempErr, par[2].tempErr, par[3].tempErr};
	int numpoints = 4;

	TGraphErrors *graph1 = new TGraphErrors(numpoints, x1, y1, ex1, ey1);
	TCanvas *c1 = new TCanvas("sca_t_txt", "sca_t_txt", 960, 570, 500, 450);

	gStyle->SetTitleY(0.96);
	graph1->SetTitle("");
	graph1->GetXaxis()->CenterTitle();
	graph1->GetXaxis()->SetTitle("Centrality");
	graph1->GetXaxis()->SetTitleOffset(1.3);
	graph1->GetYaxis()->CenterTitle();
	graph1->GetYaxis()->SetLabelSize(0.03);
	graph1->GetYaxis()->SetTitle("T_{eff}");
	graph1->GetYaxis()->SetTitleOffset(1.48);
	graph1->SetMarkerColor(4);
	graph1->SetMarkerStyle(20);
	graph1->SetMarkerSize(1.3);
	graph1->Draw("AP");
	c1->SetTickx();
	c1->SetTicky();

	TLatex *tex1 = new TLatex; tex1->SetNDC();
	tex1->SetTextSize(0.04);
	tex1->SetTextFont(42);
	tex1->DrawLatex(0.17, 0.84, "200 GeV Au+Au");
	TLatex *tex2 = new TLatex; tex2->SetNDC();
	tex2->SetTextSize(0.04);
	tex2->SetTextFont(42);
	tex2->DrawLatex(0.15, 0.15, "(c)");

	c1->SaveAs("fig1c_from_txt.eps");
}

void draw_from_parameter_txt()
{
	TsaParameter par[4];
	if (!read_parameters("tsa_parameters_one_parameter_errors.txt", par)) return;

	const int oldErrorIgnoreLevel = gErrorIgnoreLevel;
	gErrorIgnoreLevel = kFatal;
	TFile *f = TFile::Open("HEPData-ins2061074-v1-Figure_15.root", "READ");
	if (!f || f->IsZombie()) {
		gErrorIgnoreLevel = oldErrorIgnoreLevel;
		std::cerr << "Error: cannot open HEPData-ins2061074-v1-Figure_15.root\n";
		return;
	}
	TDirectoryFile *td = dynamic_cast<TDirectoryFile *>(f->Get("Figure 15"));
	if (!td) {
		gErrorIgnoreLevel = oldErrorIgnoreLevel;
		std::cerr << "Error: cannot find Figure 15\n";
		return;
	}

	draw_cenm1_from_txt(par, td);
	draw_q_from_txt(par);
	draw_t_from_txt(par);

	f->Close();
	gErrorIgnoreLevel = kFatal;
}
