#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TError.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TStyle.h"

#include <array>
#include <fstream>
#include <iostream>
#include <string>

struct FitParameterSet {
    std::string centrality;
    double nor = 0.;
    double norErrLow = 0.;
    double norErrHigh = 0.;
    double q = 0.;
    double qErrLow = 0.;
    double qErrHigh = 0.;
    double temp = 0.;
    double tempErrLow = 0.;
    double tempErrHigh = 0.;
};

Double_t nonprompt_tsallis(Double_t *x, Double_t *par)
{
    const Double_t mass = 0.;
    const Double_t b = TMath::Sqrt(TMath::Power(x[0], 2) +
                                   TMath::Power(mass, 2));
    const Double_t z1 = TMath::Power(1. + (par[1] - 1.) * b / par[2],
                                     par[1] / (1. - par[1]));
    const Double_t z2 = 2. * par[0] * b / TMath::Power(2. * TMath::Pi(), 3);
    return z1 * z2 * TMath::Power(10., 4);
}

bool read_fit_parameters(const char *path,
                         std::array<FitParameterSet, 4> &params)
{
    std::ifstream input(path);
    if (!input) {
        std::cerr << "Error: cannot open " << path << '\n';
        return false;
    }

    std::string first;
    int count = 0;
    while (input >> first) {
        if (first.empty()) continue;
        if (first[0] == '#') {
            std::string rest;
            std::getline(input, rest);
            continue;
        }
        if (count >= static_cast<int>(params.size())) {
            std::cerr << "Error: too many rows in " << path << '\n';
            return false;
        }
        auto &p = params[count];
        p.centrality = first;
        input >> p.nor >> p.norErrLow >> p.norErrHigh
              >> p.q >> p.qErrLow >> p.qErrHigh
              >> p.temp >> p.tempErrLow >> p.tempErrHigh;
        if (!input) {
            std::cerr << "Error: malformed row in " << path << '\n';
            return false;
        }
        ++count;
    }

    if (count != static_cast<int>(params.size())) {
        std::cerr << "Error: expected 4 parameter rows in " << path
                  << ", got " << count << '\n';
        return false;
    }
    return true;
}

void scale_graph(TGraphAsymmErrors *graph, double scale)
{
    for (int i = 0; i < graph->GetN(); ++i) {
        double x = 0.;
        double y = 0.;
        graph->GetPoint(i, x, y);
        graph->SetPoint(i, x, scale * y);
        graph->SetPointEYlow(i, scale * graph->GetErrorYlow(i));
        graph->SetPointEYhigh(i, scale * graph->GetErrorYhigh(i));
    }
}

void draw_entry(double y, int markerStyle, int markerColor,
                const char *label)
{
    auto *line = new TLine(0.65, y + 0.01, 0.70, y + 0.01);
    line->SetNDC();
    line->SetLineColor(kCyan - 2);
    line->SetLineWidth(2);
    line->Draw();

    auto *marker = new TMarker(0.675, y + 0.01, markerStyle);
    marker->SetNDC();
    marker->SetMarkerColor(markerColor);
    marker->SetMarkerSize(1.);
    marker->Draw("same");

    auto *text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextFont(42);
    text->DrawLatex(0.71, y, label);
}

void cenm1_nonprompt()
{
    const int oldErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    constexpr const char *rootPath = "HEPData-ins2061074-v1-Figure_16.root";
    constexpr const char *parameterPath = "qt_contours_err_results.txt";
    constexpr double xMin = 0.;
    constexpr double xMax = 12.;
    constexpr int parN = 3;

    std::array<FitParameterSet, 4> params;
    if (!read_fit_parameters(parameterPath, params)) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        return;
    }

    auto *file = TFile::Open(rootPath);
    if (!file || file->IsZombie()) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot open " << rootPath << '\n';
        return;
    }
    auto *directory = dynamic_cast<TDirectoryFile *>(file->Get("Figure 16"));
    if (!directory) {
        gErrorIgnoreLevel = oldErrorIgnoreLevel;
        std::cerr << "Error: cannot find Figure 16\n";
        return;
    }

    const std::array<const char *, 4> graphNames = {
        "Graph1D_y1", "Graph1D_y2", "Graph1D_y3", "Graph1D_y4"};
    std::array<TGraphAsymmErrors *, 4> graphs = {};
    for (std::size_t i = 0; i < graphs.size(); ++i) {
        auto *source =
            dynamic_cast<TGraphAsymmErrors *>(directory->Get(graphNames[i]));
        if (!source) {
            gErrorIgnoreLevel = oldErrorIgnoreLevel;
            std::cerr << "Error: cannot load Graph1D_y" << i + 1 << '\n';
            return;
        }
        graphs[i] = dynamic_cast<TGraphAsymmErrors *>(
            source->Clone(Form("nonprompt_graph_%zu", i + 1)));
        if (gDirectory) gDirectory->GetList()->Remove(graphs[i]);
    }

    const std::array<double, 4> scales = {30., 10., 5., 1.};
    const std::array<int, 4> markerStyles = {20, 21, 22, 23};
    const std::array<int, 4> markerColors = {4, 1, 4, 1};
    const std::array<const char *, 4> labels = {
        "0#minus20%, #times30",
        "20#minus40%, #times10",
        "40#minus60%, #times5",
        "60#minus80%, #times1",
    };

    auto *canvas = new TCanvas("c_cenm1_nonprompt", "nonprompt spectra", 960, 0,
                               550, 500);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetLogy();

    const bool oldAddToGlobalList = TF1::DefaultAddToGlobalList(false);

    gStyle->SetTitleY(0.96);
    gStyle->SetLineScalePS(2);
    gStyle->SetLineStyleString(10, "15 15");

    std::array<TF1 *, 4> funcs = {};
    for (std::size_t i = 0; i < graphs.size(); ++i) {
        scale_graph(graphs[i], scales[i]);

        graphs[i]->SetLineColor(kCyan - 2);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerStyle(markerStyles[i]);
        graphs[i]->SetMarkerColor(markerColors[i]);
        graphs[i]->SetMarkerSize(1.);

        funcs[i] = new TF1(Form("nonprompt_fit_%zu", i + 1),
                           nonprompt_tsallis, xMin, xMax, parN);
        funcs[i]->SetLineColor(2);
        funcs[i]->SetLineWidth(2);
        funcs[i]->SetParNames("Nor", "q", "T");
        funcs[i]->SetParameters(scales[i] * params[i].nor, params[i].q,
                                params[i].temp);
    }

    graphs[0]->GetXaxis()->SetLimits(0., 11.);
    graphs[0]->SetMinimum(TMath::Power(10., -9));
    graphs[0]->SetMaximum(TMath::Power(10., 3));
    graphs[0]->SetTitle("");
    graphs[0]->GetXaxis()->CenterTitle();
    graphs[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    graphs[0]->GetXaxis()->SetTitleOffset(1.3);
    graphs[0]->GetYaxis()->CenterTitle();
    graphs[0]->GetYaxis()->SetLabelSize(0.03);
    graphs[0]->GetYaxis()->SetTitle(
        "1/N 1/(2#pip_{T}) d^{2}N/(dydp_{T}) ((GeV/c)^{-2})");
    graphs[0]->GetYaxis()->SetTitleOffset(1.3);

    graphs[0]->Draw("AP");
    funcs[0]->Draw("same");
    for (std::size_t i = 1; i < graphs.size(); ++i) {
        graphs[i]->Draw("P same");
        funcs[i]->Draw("same");
    }

    auto *fitLine = new TLine(0.20, 0.78, 0.25, 0.78);
    fitLine->SetNDC();
    fitLine->SetLineStyle(1);
    fitLine->SetLineColor(2);
    fitLine->SetLineWidth(2);
    fitLine->Draw("same");

    auto *fitText = new TLatex;
    fitText->SetNDC();
    fitText->SetTextSize(0.03);
    fitText->SetTextFont(42);
    fitText->DrawLatex(0.27, 0.77, "Tsallis");

    auto *title = new TLatex;
    title->SetNDC();
    title->SetTextSize(0.03);
    title->SetTextFont(132);
    title->DrawLatex(0.2, 0.84, "nonprompt-photon 200 GeV Au+Au, PHENIX");

    draw_entry(0.84, markerStyles[0], markerColors[0], labels[0]);
    draw_entry(0.79, markerStyles[1], markerColors[1], labels[1]);
    draw_entry(0.74, markerStyles[2], markerColors[2], labels[2]);
    draw_entry(0.69, markerStyles[3], markerColors[3], labels[3]);

    auto *panel = new TLatex;
    panel->SetNDC();
    panel->SetTextSize(0.04);
    panel->SetTextFont(42);
    panel->DrawLatex(0.15, 0.15, "(a)");

    canvas->SaveAs("cenm1_nonprompt.png");
    canvas->SaveAs("cenm1_nonprompt.pdf");

    TF1::DefaultAddToGlobalList(oldAddToGlobalList);
}
