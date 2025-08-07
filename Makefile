.PHONY: all qc_reads scan_sequences map_plasmid clean

all: qc_reads scan_sequences map_plasmid

qc_reads:
	python scripts/fastp.py -c configs.yaml

scan_sequences:
	python scripts/scan_sequences.py -c configs.yaml

map_plasmid:
	python scripts/map_to_plasmid.py -c configs.yaml

clean:
	rm -rf results/*
	rm -rf inputs/qc_reads/*
	rm -rf inputs/paired/*
	find . -name "*.log" -delete
	find . -name "*.tmp" -delete
