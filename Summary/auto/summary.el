(TeX-add-style-hook "summary"
 (lambda ()
    (LaTeX-add-labels
     "figOneFiber"
     "figTwoFibers1"
     "figTwoFibers2"
     "figPhantomDiffus"
     "figPhantomWeights")
    (TeX-run-style-hooks
     "float"
     "subcaption"
     "caption"
     "graphicx"
     "latex2e"
     "art10"
     "article")))

