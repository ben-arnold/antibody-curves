# antibody-curves

This repository includes scripts required to conduct all of the analyses for the article entitled:

_Measuring changes in transmission of neglected tropical diseases, malaria, and enteric pathogens from quantitative antibody levels_

Published in _PLoS Neglected Tropical Diseases_ in 2017:
[http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005616](http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005616)

### Notes

The sub-directories are organized by datasource (enterics, garki, mauke). Within each directory, scripts will be named according to the figure or table that they generate. The analyses rely a lot on the [`tmleAb`](http://www.github.com/ben-arnold/tmleAb) R package. In fact, the package includes all of the datasets as well! We have also archived the scripts (and datasets in .csv format) through the Open Science Framework ([https://osf.io/8tqu4](https://osf.io/8tqu4/)).  

Note that this repo includes a script directory for Miton, Haiti, which was another malaria example that we excluded from the final paper submitted to PLOS only because of space constraints (the Miton data are included in the `tmleAb` package as [`miton_malaria`](http://ben-arnold.github.io/tmleAb/miton_malaria.html)).