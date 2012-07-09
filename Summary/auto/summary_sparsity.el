(TeX-add-style-hook "summary_sparsity"
 (lambda ()
    (LaTeX-add-bibliographies
     "ref")
    (LaTeX-add-labels
     "fig:cc"
     "fig:diffusivities"
     "fig:sparsity_roi"
     "fig:sparsity_plot"
     "fig:roi")
    (TeX-run-style-hooks
     "chngpage"
     "subcaption"
     "caption"
     "graphicx"
     "latex2e"
     "art10"
     "article")))

