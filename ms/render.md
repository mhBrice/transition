

## Main text

pandoc ms/Brice_ms_transition.md -o ms/Brice_ms_transition.pdf  --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

pandoc ms/Brice_ms_transition.md -f markdown -t latex -s -o ms/Brice_ms_transition.tex --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

pandoc ms/Brice_ms_transition.md -o ms/Brice_ms_transition.docx  --bibliography=../references.bib --csl ms/GCB.csl

## Supplementary material

pandoc ms/Brice_SI_transition.md -o ms/Brice_SI_transition.pdf --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex
