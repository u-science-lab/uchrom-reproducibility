# The manuscripts' figure files and generated tables, then copied into the paper checkout.  Draw a figure with its
# make_*.py first (each says what it reads); `make sync` copies what the paper tracks (sync_to_paper.py).
.PHONY: full tables sync
FULL = manuscripts/full

# The full manuscript's figures as its LaTeX includes them (formerly `make figures` in the paper's manuscripts/full).
full:
	cp $(FULL)/figures/fig1/fig1_final.pdf $(FULL)/figures/fig1.pdf
	cp $(FULL)/figures/fig2/fig2.pdf $(FULL)/figures/fig2.pdf
	cp $(FULL)/figures/fig2/figS_fig2.pdf $(FULL)/figures/figS1.pdf
	cp $(FULL)/figures/fig6/fig6.pdf $(FULL)/figures/fig6.pdf

# Supplementary Table S1's bodies (full) and the atlas table (nmeth-bc), from their results.
tables:
	cd $(FULL) && python make_tableS1.py figures/fig2/tableS_capabilities.md supp
	cd manuscripts/nmeth-bc && python make_supp_tables.py

sync:
	python sync_to_paper.py
