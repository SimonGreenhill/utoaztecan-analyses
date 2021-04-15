
TREES = ../phylogeny/ua-covarion-relaxed.mcct.trees
VARIABLES = EA005 EA028 EA029 EA042
DATAFILES = $(patsubst %, %.tsv, $(VARIABLES))
RESULTFILES = $(patsubst %, %.result.csv, $(VARIABLES))
PDFFILES = $(patsubst %, %.pdf, $(VARIABLES))

all: $(RESULTFILES) $(PDFFILES)

%.result.csv: %.tsv
	./get_rates.R $(TREES) $< $@

%.pdf: %.result.csv
	./plot.R $< $@ > $(patsubst %.result.csv, %.result.txt, $<)

clean:
	rm -rf *.result.csv *.pdf
