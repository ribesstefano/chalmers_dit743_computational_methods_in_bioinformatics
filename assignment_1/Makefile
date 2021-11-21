CXX := gcc
CFLAGS := -std=c99

.PHONY: all

all: global_alignment.exe levenshtein.exe local_alignment.exe

global_alignment.exe: global_alignment.c
	$(CXX) $(CFLAGS) $< -o $@

levenshtein.exe: levenshtein.c
	$(CXX) $(CFLAGS) $< -o $@

local_alignment.exe: local_alignment.c
	$(CXX) $(CFLAGS) $< -o $@


run: global_alignment.exe levenshtein.exe local_alignment.exe
	@echo "# ========================================================"
	@echo "# Running Global Alignment"
	@echo "# ========================================================"
	./global_alignment.exe
	@echo "# ========================================================"
	@echo "# Running Levenshtein"
	@echo "# ========================================================"
	./levenshtein.exe
	@echo "# ========================================================"
	@echo "# Running Local Alignment"
	@echo "# ========================================================"
	./local_alignment.exe
	@echo "# ========================================================"
	@echo "# Running Optimal Paths Search"
	@echo "# ========================================================"
	./global_alignment.exe ATTA ATTTTA

clean:
	rm -rf *.exe *.o
