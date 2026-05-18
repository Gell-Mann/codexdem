#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct Options {
  int events = 1000;
  int seed = 12345;
  std::string config = "jpsi_200gev.cmnd";
  std::string hepmc = "jpsi_200gev.hepmc";
  std::string root = "jpsi_pt_spectrum.root";
  std::string pdf = "jpsi_pt_spectrum.pdf";
  std::string dat = "jpsi_pt_spectrum.dat";
};

int toInt(const std::string &value, const std::string &name) {
  char *end = nullptr;
  const long parsed = std::strtol(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    throw std::runtime_error("Invalid integer for " + name + ": " + value);
  }
  return static_cast<int>(parsed);
}

Options parseArgs(int argc, char **argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto requireValue = [&](const std::string &name) -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error("Missing value after " + name);
      }
      return argv[++i];
    };

    if (arg == "--events" || arg == "-n") {
      opt.events = toInt(requireValue(arg), arg);
    } else if (arg == "--seed") {
      opt.seed = toInt(requireValue(arg), arg);
    } else if (arg == "--config") {
      opt.config = requireValue(arg);
    } else if (arg == "--hepmc") {
      opt.hepmc = requireValue(arg);
    } else if (arg == "--root") {
      opt.root = requireValue(arg);
    } else if (arg == "--pdf") {
      opt.pdf = requireValue(arg);
    } else if (arg == "--dat") {
      opt.dat = requireValue(arg);
    } else if (arg == "--help" || arg == "-h") {
      std::cout
          << "Usage: " << argv[0] << " [options]\n"
          << "  -n, --events N       Number of generated events [1000]\n"
          << "      --seed N         Pythia random seed [12345]\n"
          << "      --config FILE    Extra Pythia command file [jpsi_200gev.cmnd]\n"
          << "      --hepmc FILE     HepMC3 ASCII output [jpsi_200gev.hepmc]\n"
          << "      --root FILE      ROOT histogram output [jpsi_pt_spectrum.root]\n"
          << "      --pdf FILE       Spectrum plot output [jpsi_pt_spectrum.pdf]\n"
          << "      --dat FILE       Binned text output [jpsi_pt_spectrum.dat]\n";
      std::exit(EXIT_SUCCESS);
    } else {
      throw std::runtime_error("Unknown option: " + arg);
    }
  }

  if (opt.events <= 0) {
    throw std::runtime_error("Number of events must be positive");
  }
  if (opt.seed <= 0 || opt.seed > 900000000) {
    throw std::runtime_error("Pythia seed must be in the range 1..900000000");
  }
  return opt;
}

void readSetting(Pythia8::Pythia &pythia, const std::string &setting) {
  if (!pythia.readString(setting)) {
    throw std::runtime_error("Pythia rejected setting: " + setting);
  }
}

bool fileExists(const std::string &path) {
  std::ifstream in(path);
  return in.good();
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Options opt = parseArgs(argc, argv);

    Pythia8::Pythia pythia;
    readSetting(pythia, "Beams:idA = 2212");
    readSetting(pythia, "Beams:idB = 2212");
    readSetting(pythia, "Beams:eCM = 200.");
    readSetting(pythia, "Random:setSeed = on");
    readSetting(pythia, "Random:seed = " + std::to_string(opt.seed));
    readSetting(pythia, "Next:numberShowInfo = 0");
    readSetting(pythia, "Next:numberShowProcess = 0");
    readSetting(pythia, "Next:numberShowEvent = 0");

    readSetting(pythia, "Charmonium:all = on");
    readSetting(pythia, "443:mayDecay = off");

    if (!opt.config.empty() && fileExists(opt.config)) {
      if (!pythia.readFile(opt.config)) {
        throw std::runtime_error("Failed to read Pythia command file: " +
                                 opt.config);
      }
    }

    if (!pythia.init()) {
      throw std::runtime_error("Pythia initialization failed");
    }

    Pythia8::Pythia8ToHepMC toHepMC(opt.hepmc);

    TH1D hPt("hJpsiPt", "J/#psi p_{T} in pp at #sqrt{s}=200 GeV;p_{T} [GeV/c];Counts", 100,
             0.0, 20.0);
    hPt.Sumw2();

    long acceptedEvents = 0;
    long jpsiCount = 0;

    for (int iEvent = 0; iEvent < opt.events; ++iEvent) {
      if (!pythia.next()) {
        continue;
      }
      ++acceptedEvents;

      for (int i = 0; i < pythia.event.size(); ++i) {
        const Pythia8::Particle &p = pythia.event[i];
        if (p.id() == 443 && p.isFinal()) {
          hPt.Fill(p.pT());
          ++jpsiCount;
        }
      }

      if (!toHepMC.writeNextEvent(pythia)) {
        throw std::runtime_error("Failed to write HepMC event");
      }
    }

    TFile outFile(opt.root.c_str(), "RECREATE");
    hPt.Write();
    outFile.Close();

    gStyle->SetOptStat(0);
    TCanvas canvas("cJpsiPt", "J/psi pT spectrum", 800, 600);
    hPt.SetLineColor(kBlue + 1);
    hPt.SetMarkerColor(kBlue + 1);
    hPt.SetMarkerStyle(20);
    hPt.Draw("E1");
    canvas.SaveAs(opt.pdf.c_str());

    std::ofstream datOut(opt.dat);
    datOut << "# bin_low_GeV bin_high_GeV bin_center_GeV counts stat_error\n";
    for (int bin = 1; bin <= hPt.GetNbinsX(); ++bin) {
      datOut << hPt.GetBinLowEdge(bin) << ' '
             << hPt.GetBinLowEdge(bin + 1) << ' '
             << hPt.GetBinCenter(bin) << ' '
             << hPt.GetBinContent(bin) << ' '
             << hPt.GetBinError(bin) << '\n';
    }

    pythia.stat();
    std::cout << "Generated events requested: " << opt.events << '\n'
              << "Accepted Pythia events: " << acceptedEvents << '\n'
              << "Final-state J/psi count: " << jpsiCount << '\n'
              << "Wrote HepMC: " << opt.hepmc << '\n'
              << "Wrote ROOT histogram: " << opt.root << '\n'
              << "Wrote PDF spectrum: " << opt.pdf << '\n'
              << "Wrote binned table: " << opt.dat << '\n';
  } catch (const std::exception &err) {
    std::cerr << "Error: " << err.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
