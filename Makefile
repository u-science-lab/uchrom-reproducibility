# Generated tables, then the copy into the paper checkout.  Draw a figure with its script first (each folder's README
# says what it reads); `make sync` copies what the paper includes (sync_to_paper.py).
.PHONY: tables sync

tables:
	python tables/supp_table2/make_supp_tables.py

sync:
	python sync_to_paper.py
