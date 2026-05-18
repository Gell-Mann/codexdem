CXX ?= g++

LOCAL_PYTHIA8_CONFIG = /home/weyl/pythia8/pythia8307/bin/pythia8-config
PYTHIA8_CONFIG ?= $(if $(wildcard $(LOCAL_PYTHIA8_CONFIG)),$(LOCAL_PYTHIA8_CONFIG),pythia8-config)
HEPMC3_CONFIG ?= HepMC3-config
ROOT_CONFIG ?= root-config

CXXFLAGS_BASE = -O2 -std=c++17 -Wall -Wextra -pedantic
CXXFLAGS_DEPS = $(shell $(PYTHIA8_CONFIG) --cxxflags) \
                $(shell $(HEPMC3_CONFIG) --cflags) \
                $(shell $(ROOT_CONFIG) --cflags)

LDLIBS_DEPS = $(shell $(PYTHIA8_CONFIG) --libs) \
              $(shell $(HEPMC3_CONFIG) --libs) \
              $(shell $(ROOT_CONFIG) --libs)

TARGET = simulate_jpsi_200gev
SRC = simulate_jpsi_200gev.cc

.PHONY: all check-deps clean run

all: $(TARGET)

check-deps:
	@command -v $(PYTHIA8_CONFIG) >/dev/null || { echo "Missing $(PYTHIA8_CONFIG). Install Pythia8 or set PYTHIA8_CONFIG=/path/to/pythia8-config."; exit 1; }
	@command -v $(HEPMC3_CONFIG) >/dev/null || { echo "Missing $(HEPMC3_CONFIG). Install HepMC3 or set HEPMC3_CONFIG=/path/to/HepMC3-config."; exit 1; }
	@command -v $(ROOT_CONFIG) >/dev/null || { echo "Missing $(ROOT_CONFIG). Install ROOT or set ROOT_CONFIG=/path/to/root-config."; exit 1; }

$(TARGET): check-deps $(SRC)
	$(CXX) $(CXXFLAGS_BASE) $(CXXFLAGS_DEPS) $(SRC) -o $@ $(LDLIBS_DEPS)

run: $(TARGET)
	./$(TARGET) --events 1000

clean:
	rm -f $(TARGET) *.o *.root *.pdf *.dat *.hepmc
