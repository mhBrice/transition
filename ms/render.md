

## Main text

pandoc ms/Brice_ms_transition.md -o ms/Brice_ms_transition.pdf  --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex


## Supplementary material

pandoc ms/revision/SI_revised.md -o ms/revision/Brice_SI_revised.pdf --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex
