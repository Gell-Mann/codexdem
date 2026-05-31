#!/usr/bin/env python3
"""Fit the extracted pi0 pT spectrum with a Tsallis blast-wave function.

This is a Python/PyROOT translation of the TBW function in
/home/weyl/rodemo/TBW/TBW.C, adapted to read the CSV spectrum produced from
UrQMD f13 output.

The script is intentionally not executed here. Typical usage:

    python3 fit_pi0_tbw_from_csv.py \
      --spectrum /home/weyl/codexdem/pi0_f13_pt_spectrum.csv \
      --events 10

The y-value is built from the binned counts as
    (1 / N_event) * dN / (2*pi*pT*dpT)
with Poisson errors. This is not dN/dy-normalized because the current UrQMD
CSV was not rapidity-selected.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import ROOT


PI0_MASS_GEV = 0.139
TBW_RADIUS_FM = 22.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit pi0 pT spectrum CSV with TBW.")
    parser.add_argument(
        "--spectrum",
        default="/home/weyl/codexdem/pi0_f13_pt_spectrum.csv",
        help="CSV with columns pt_low_GeV_c, pt_high_GeV_c, pt_center_GeV_c, count.",
    )
    parser.add_argument(
        "--events",
        type=float,
        default=10.0,
        help="Number of UrQMD events used to normalize the spectrum.",
    )
    parser.add_argument("--fit-min", type=float, default=0.0, help="Lower pT fit limit.")
    parser.add_argument("--fit-max", type=float, default=4.0, help="Upper pT fit limit.")
    parser.add_argument("--mass", type=float, default=PI0_MASS_GEV, help="Particle mass in GeV.")
    parser.add_argument("--radius", type=float, default=TBW_RADIUS_FM, help="TBW radius parameter.")
    parser.add_argument(
        "--fix-beta",
        type=float,
        default=0.37,
        help="Fix TBW beta parameter to this value. Default: 0.37.",
    )
    parser.add_argument("--nr", type=int, default=80, help="Number of radial grid points.")
    parser.add_argument("--nphi", type=int, default=72, help="Number of azimuthal grid points.")
    parser.add_argument("--out-root", default="/home/weyl/codexdem/pi0_tbw_fit_pt0_4_beta037.root")
    parser.add_argument("--out-pdf", default="/home/weyl/codexdem/pi0_tbw_fit_pt0_4_beta037.pdf")
    parser.add_argument("--out-txt", default="/home/weyl/codexdem/pi0_tbw_fit_pt0_4_beta037_summary.txt")
    return parser.parse_args()


def read_spectrum(path: str, events: float, fit_min: float, fit_max: float):
    points = []
    with Path(path).open("r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin)
        for row in reader:
            try:
                pt_low = float(row["pt_low_GeV_c"])
                pt_high = float(row["pt_high_GeV_c"])
                pt = float(row["pt_center_GeV_c"])
                count = float(row["count"])
            except ValueError:
                continue

            if pt <= 0.0 or count <= 0.0:
                continue
            if pt < fit_min or pt > fit_max:
                continue

            width = pt_high - pt_low
            yield_value = count / (events * 2.0 * math.pi * pt * width)
            yield_error = math.sqrt(count) / (events * 2.0 * math.pi * pt * width)
            x_error = 0.5 * width
            points.append((pt, yield_value, x_error, yield_error, count))
    return points


class TBWFunction:
    """Callable PyROOT TF1 body matching the TBW.C formula.

    Parameters:
        par[0] T0   kinetic temperature-like parameter, GeV
        par[1] n0   radial profile exponent
        par[2] beta maximum transverse flow velocity
        par[3] q    Tsallis q parameter
        par[4] Nor  normalization
    """

    def __init__(self, mass: float, radius: float, nr: int, nphi: int):
        self.mass = mass
        self.radius = radius
        self.nr = nr
        self.nphi = nphi
        self.dr = radius / nr
        self.dphi = 2.0 * math.pi / nphi
        self.r_midpoints = [(i + 0.5) * self.dr for i in range(nr)]
        self.phi_midpoints = [-math.pi + (j + 0.5) * self.dphi for j in range(nphi)]
        self.cos_phi = [math.cos(phi) for phi in self.phi_midpoints]

    def __call__(self, x, par):
        pt = float(x[0])
        temp = float(par[0])
        n0 = float(par[1])
        beta = float(par[2])
        q = float(par[3])
        norm = float(par[4])

        if temp <= 0.0 or n0 <= 0.0 or beta <= 0.0 or beta >= 1.0 or q <= 1.0:
            return 0.0

        mt = math.sqrt(pt * pt + self.mass * self.mass)
        q_minus_one = q - 1.0
        exponent = -1.0 / q_minus_one
        integral = 0.0

        for r in self.r_midpoints:
            br = beta * (r / self.radius) ** n0
            if br >= 1.0:
                return 0.0
            rho = math.atanh(br)
            sinh_rho = math.sinh(rho)
            cosh_rho = math.cosh(rho)
            radial_sum = 0.0
            for cos_phi in self.cos_phi:
                energy_shift = mt * cosh_rho - pt * sinh_rho * cos_phi
                base = 1.0 + q_minus_one * energy_shift / temp
                if base <= 0.0:
                    continue
                radial_sum += base**exponent
            integral += r * radial_sum * self.dphi * self.dr

        return mt * integral * norm


def make_graph(points):
    graph = ROOT.TGraphErrors(len(points))
    for i, (pt, y, ex, ey, _count) in enumerate(points):
        graph.SetPoint(i, pt, y)
        graph.SetPointError(i, ex, ey)
    graph.SetName("pi0_pt_spectrum")
    graph.SetTitle("UrQMD #pi^{0} p_{T} spectrum; p_{T} (GeV/c); (1/N_{evt}) dN/(2#pi p_{T} dp_{T})")
    return graph


def main() -> int:
    args = parse_args()
    points = read_spectrum(args.spectrum, args.events, args.fit_min, args.fit_max)
    if len(points) < 5:
        raise SystemExit("Not enough nonzero pT bins in the requested fit range.")

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptFit(1111)

    graph = make_graph(points)
    tbw_callable = TBWFunction(args.mass, args.radius, args.nr, args.nphi)
    fit_func = ROOT.TF1("tbw_pi0", tbw_callable, args.fit_min, args.fit_max, 5)
    fit_func.SetParNames("T0", "n0", "beta", "q", "Nor")

    fit_func.SetParameter(0, 0.089)
    fit_func.SetParameter(1, 2.0)
    fit_func.SetParameter(2, args.fix_beta if args.fix_beta is not None else 0.83)
    fit_func.SetParameter(3, 1.04)
    fit_func.SetParameter(4, 1.0e-3)

    fit_func.SetParLimits(0, 0.001, 1.0)
    fit_func.FixParameter(1, 2.0)
    if args.fix_beta is None:
        fit_func.SetParLimits(2, 0.01, 0.999)
    else:
        fit_func.FixParameter(2, args.fix_beta)
    fit_func.SetParLimits(3, 1.0001, 1.2)
    fit_func.SetParLimits(4, 1.0e-12, 1.0e6)

    fit_result = graph.Fit(fit_func, "S R")

    canvas = ROOT.TCanvas("c_pi0_tbw", "pi0 TBW fit", 800, 650)
    canvas.SetLogy(True)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.Draw("AP")
    fit_func.SetLineColor(ROOT.kRed + 1)
    fit_func.Draw("same")
    canvas.SaveAs(args.out_pdf)

    out_file = ROOT.TFile(args.out_root, "RECREATE")
    graph.Write()
    fit_func.Write()
    canvas.Write()
    out_file.Close()

    with Path(args.out_txt).open("w", encoding="utf-8") as fout:
        fout.write(f"spectrum: {args.spectrum}\n")
        fout.write(f"events: {args.events:g}\n")
        fout.write(f"fit range GeV/c: {args.fit_min:g} to {args.fit_max:g}\n")
        fout.write(f"mass GeV: {args.mass:g}\n")
        fout.write(f"radius fm: {args.radius:g}\n")
        fout.write(f"fixed beta: {args.fix_beta if args.fix_beta is not None else 'None'}\n")
        fout.write(f"integration grid: nr={args.nr} nphi={args.nphi}\n")
        fout.write(f"fit status: {int(fit_result.Status())}\n")
        fout.write(f"chi2: {fit_func.GetChisquare():.10g}\n")
        fout.write(f"ndf: {fit_func.GetNDF()}\n")
        for ipar in range(5):
            fout.write(
                f"{fit_func.GetParName(ipar)}: "
                f"{fit_func.GetParameter(ipar):.10g} +/- {fit_func.GetParError(ipar):.10g}\n"
            )
        fout.write(f"root output: {args.out_root}\n")
        fout.write(f"pdf output: {args.out_pdf}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
