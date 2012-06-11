(TeX-add-style-hook "summary"
 (lambda ()
    (LaTeX-add-labels
     "figOneFiber")
    (TeX-run-style-hooks
     "subcaption"
     "caption"
     "graphicx"
     "latex2e"
     "art10"
     "article")))

