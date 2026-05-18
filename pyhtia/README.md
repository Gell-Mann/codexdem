# Pythia8 J/psi production at 200 GeV pp

This directory contains a small Pythia8 driver for generating `pp` events at
`sqrt(s) = 200 GeV`, selecting final-state `J/psi` (`PDG 443`) particles, and
writing:

- `jpsi_200gev.hepmc`: generated events in HepMC3 ASCII format
- `jpsi_pt_spectrum.root`: ROOT histogram with the `J/psi` transverse momentum
- `jpsi_pt_spectrum.pdf`: plotted `pT` spectrum
- `jpsi_pt_spectrum.dat`: binned histogram table

The generator settings are intentionally plain and editable. They are suitable
for producing a simulation sample, not for claiming a tuned physical result.

## Build

Pythia8, HepMC3, and ROOT must be available in the environment:

```bash
make
```

The Makefile uses `pythia8-config`, `HepMC3-config`, and `root-config`.

## Run

```bash
./simulate_jpsi_200gev --events 1000
```

Useful options:

```bash
./simulate_jpsi_200gev \
  --events 100000 \
  --seed 12345 \
  --hepmc jpsi_200gev.hepmc \
  --root jpsi_pt_spectrum.root \
  --pdf jpsi_pt_spectrum.pdf \
  --dat jpsi_pt_spectrum.dat
```

Extra Pythia settings can be placed in `jpsi_200gev.cmnd` and loaded with:

```bash
./simulate_jpsi_200gev --config jpsi_200gev.cmnd
```
