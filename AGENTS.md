# Project AGENTS.md

## Project goal
This repository is for high-energy / nuclear physics analysis involving:
- transverse momentum spectra fitting with ROOT;
- event generation with Pythia and Herwig;
- detector simulation with Geant4;
- jet-related ML studies with ROOT TMVA and Python;
- neutron-star nuclear symmetry-energy analysis.

## Repository layout
- `src/`: C++ and Python source code.
- `macros/`: ROOT macros.
- `scripts/`: reproducible analysis scripts.
- `config/pythia/`: Pythia command files.
- `config/herwig/`: Herwig input files.
- `config/geant4/`: Geant4 macro/configuration files.
- `docs/`: project conventions and physics notes.
- `tests/`: small validation tests.
- `results/`: generated figures, tables, and logs.
- `data/`: large external data. Do not inspect in full.

## Default commands
- Build C++ code: `cmake -S . -B build && cmake --build build -j`
- Run tests: `ctest --test-dir build --output-on-failure`
- Format Python: `ruff format .`
- Check Python: `ruff check .`
- Run ROOT batch macro: `root -l -b -q macros/example.C`
- Prefer small validation jobs before full production jobs.

## Physics rules
- Always track units.
- Always track random seeds and software versions.
- Do not compare spectra, cross sections, efficiencies, or yields unless normalization conventions are explicit.
- For fits, always report fit function, fit range, bin width treatment, parameter uncertainties, chi2/ndf or likelihood metric, and residual behavior.
- For MC generation, always record collision system, sqrt(s), tune, PDF, process switches, cuts, event count, seed, and cross-section normalization.
- For Geant4, always record geometry version, material definition, physics list, production cuts, sensitive detector mapping, and random seed.
- For ML, always check train/test split, leakage, feature definitions, preprocessing, class imbalance, overtraining, and systematic robustness.
- For neutron-star EOS/symmetry energy, check density units, pressure monotonicity, causality, beta equilibrium, charge neutrality, and TOV reproducibility where relevant.

## Token-saving rules
- Do not open large ROOT files directly. Write or use scripts that print compact metadata.
- Do not dump full logs. Use `tail -n 100`, grep relevant errors, or produce compact summaries.
- Do not paste long papers or manuals into context. Use references or MCP lookup only when needed.
- Prefer modifying existing files over generating extensive new boilerplate.

## Definition of done
A task is done only when:
1. the requested code or analysis change is completed;
2. relevant minimal checks have been run or clearly explained if not run;
3. assumptions and remaining physical/software risks are listed;
4. output files, figures, or tables are named explicitly if produced.
