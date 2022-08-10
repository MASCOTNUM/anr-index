# anr-index

This repository contains software developped in the
[ANR project INDEX](https://sdb3.i3s.unice.fr/anrindex/).

![ANR INDEX logo](anr-index.png)


## Content

* Incremental construction of designs that greedily maximize the criterion I_B,q in the paper "Incremental space-filling design based on coverings and spacings: improving upon low discrepancy sequences", by Nogales-Gómez, Pronzato and Rendas, submitted to the Journal of Statistical Theory and Practice (2021):
  * [rungreedy_cdf.m](mfiles/rungreedy_cdf.m)
  * [computei_a.m](mfiles/computei_a.m)
  * [greedy_cdf.m](mfiles/greedy_cdf.m)
  * [lazy_greedy_cdf.m](mfiles/lazy_greedy_cdf.m)
* Construction of coffee-house designs and edgephobe coffee-house designs in the same paper:
  * [coffee_house.m](mfiles/coffee_house.m)
* Calculation of the covering radius and covering alpha-quantile of a given design (finite approximation):
  * [covering_radius.m](mfiles/covering_radius.m)
  * [covering_quantile.m](mfiles/covering_quantile.m)
* Graphic output:
  * [numlabelplot.m](mfiles/numlabelplot.m)


## Authors

[Luc Pronzato](https://www.i3s.unice.fr/lpronzato/) &
[Maria João Rendas](https://www.i3s.unice.fr/~rendas/rendas/Home.html)


## Licence

Copyright 2021 CNRS

This software is distributed under the "simplified BSD" (BSD-2-Clause) license,
see [LICENSE](LICENSE) for details.
