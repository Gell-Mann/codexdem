# Rivet comparison for `jpsi_200gev.hepmc`

Input HepMC file:

- `../jpsi_200gev.hepmc`

Rivet setup found under:

- `/home/weyl/rivet/local/rivetenv.sh`
- `/home/weyl/rivet/local/bin/rivet`
- Rivet version: 3.1.7

## Analysis used

The generated event file is `pp` at `sqrt(s) = 200 GeV` with charmonium
production enabled. The available Rivet installation does not contain a
validated `pp 200 GeV J/psi pT` reference analysis.

The energy-matched analysis used for the produced comparison is:

- `STAR_2008_S7869363`
- Observable: charged multiplicity and identified charged hadron pT spectra
- Experiment: STAR at RHIC
- Beams: `p+ p+`, `sqrt(s) = 200 GeV`
- Rivet status: UNVALIDATED

This is a comparison of the full generated event sample to STAR charged-hadron
reference data. It is not a direct experimental validation of the generated
`J/psi` pT spectrum.

The J/psi-specific Rivet analysis `ATLAS_2011_S9035664` was checked, but Rivet
rejected it for this input because that analysis expects `pp` beams at
`sqrt(s) = 7 TeV`, not 200 GeV.

## Commands run

```bash
rivet -a STAR_2008_S7869363 \
  -o rivet_comparison/star_2008_vs_jpsi_pythia.yoda \
  jpsi_200gev.hepmc

rivet-mkhtml \
  rivet_comparison/star_2008_vs_jpsi_pythia.yoda \
  --outputdir rivet_comparison/html
```

## Output files

Main Rivet output:

- `star_2008_vs_jpsi_pythia.yoda`

Browsable comparison page:

- `html/index.html`
- `html/STAR_2008_S7869363/index.html`

Per-observable comparison data and plots:

- `html/STAR_2008_S7869363/d01-x01-y01.dat`
- `html/STAR_2008_S7869363/d01-x01-y01.png`
- `html/STAR_2008_S7869363/d01-x01-y01.pdf`
- `html/STAR_2008_S7869363/d02-x01-y01.dat`
- `html/STAR_2008_S7869363/d02-x01-y01.png`
- `html/STAR_2008_S7869363/d02-x01-y01.pdf`
- `html/STAR_2008_S7869363/d02-x01-y02.dat`
- `html/STAR_2008_S7869363/d02-x01-y02.png`
- `html/STAR_2008_S7869363/d02-x01-y02.pdf`
- `html/STAR_2008_S7869363/d02-x01-y03.dat`
- `html/STAR_2008_S7869363/d02-x01-y03.png`
- `html/STAR_2008_S7869363/d02-x01-y03.pdf`
- `html/STAR_2008_S7869363/d02-x01-y04.dat`
- `html/STAR_2008_S7869363/d02-x01-y04.png`
- `html/STAR_2008_S7869363/d02-x01-y04.pdf`
- `html/STAR_2008_S7869363/d02-x01-y05.dat`
- `html/STAR_2008_S7869363/d02-x01-y05.png`
- `html/STAR_2008_S7869363/d02-x01-y05.pdf`
- `html/STAR_2008_S7869363/d02-x01-y06.dat`
- `html/STAR_2008_S7869363/d02-x01-y06.png`
- `html/STAR_2008_S7869363/d02-x01-y06.pdf`

The `.dat` files contain both the Rivet reference data and the generated-sample
histogram. The generated plots include legends with the reference data and the
`star_2008_vs_jpsi_pythia` sample.

