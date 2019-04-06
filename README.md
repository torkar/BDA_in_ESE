# BDA_in_ESE
Reproducibility package for the chapter in Contemporary Empirical Methods in Software Engineering

There are two files in this repository: `missing_data_example.R` and `toy_example.R`.

`toy_example.R` contains the toy example which we use to allow the reader some insight into Bayesian data analysis (Sect. 2 of the chapter). Data for this example is available through the `Rethinking`package in `R`. In the example we use both the `rethinking` way of declaring models, and the `brms` way (which is based on lme4 syntax).

`missing_data_example.R` contains the actual case study (Sect. 3) where we also deal with missing data. In this case the data is not available to the reader, but instead one must buy a license for it at [ISBSG](https://www.isbsg.org).
