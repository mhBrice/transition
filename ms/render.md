
## Title page

pandoc ms/title_page.md -o ms/reviews/title_page.pdf --pdf-engine=xelatex

## Main text

pandoc ms/Brice_ms_transition.md -o ms/submission/Brice_ms_transition.pdf --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

pandoc ms/Brice_ms_transition.md -f markdown -t latex -s -o ms/Brice_ms_transition.tex --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

pandoc ms/Brice_ms_transition.md -o ms/submission/Brice_ms_transition.docx  --bibliography=../references.bib --csl ms/GCB.csl

## Supplementary material

pandoc ms/Brice_SI_transition.md -o ms/submission/Brice_SI_transition.pdf --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

## Submission

pandoc ms/GCB_submission_checklist.md -o ms/submission/GCB_submission_checklist.docx  --bibliography=../references.bib --csl ms/GCB.csl

## Revisions

pandoc ms/reviews/response2reviews.md -o ms/reviews/response2reviews.pdf  --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

pandoc ms/Brice_SI_transition.md -o ms/reviews/Brice_SI_transition.pdf --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex

pandoc ms/Brice_ms_transition.md -o ms/reviews/Brice_ms_transition.pdf  --bibliography=../references.bib --csl ms/GCB.csl --pdf-engine=xelatex
