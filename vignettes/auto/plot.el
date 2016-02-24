(TeX-add-style-hook
 "plot"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("mathpazo" "sc") ("placeins" "section")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "graphicx"
    "color"
    "framed"
    "alltt"
    "mathtools"
    "mathpazo"
    "geometry"
    "morefloats"
    "placeins"
    "longtable"
    "booktabs")))

