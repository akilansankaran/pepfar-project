# tables/

Every output the analysis scripts write. One folder per study (Crown = economic outcomes, Kates = mortality) and type (DiD or event studies). Never edit these by hand: rerun the script listed and they regenerate.

```
tables/
├── crown_did/                      Crown et al. (2023) difference-in-differences
│   ├── pooled_ols/                 Table A8-style DiD, unadjusted + adjusted   ← analysis/replication/did_crown.R  (Rishad)
│   ├── adjusted/                   adjusted pooled OLS + TWFE DiD              ← analysis/replication/did_econ.R   (Anthony)
│   └── transformations/            DiD + parallel-trends tests with raw / winsorized / log(x+1) / IHS outcomes ← analysis/replication/did_transformations.R (Anthony)
├── crown_event_studies/            Crown economic outcomes, year-by-year
│   ├── adjusted/                   adjusted event studies: tex/ figures/ csv/  ← analysis/parallel_trends/event_study_econ_adjusted.R (Anthony)
│   └── robustness/                 unadjusted, linear trends, Sun–Abraham, etc. ← analysis/parallel_trends/event_study_economic.R      (Anthony)
│                                   tex/linear_trends_*.tex                       ← analysis/parallel_trends/linear_trends_table.R     (Anthony)
├── kates_event_studies/            Kates et al. (2026) mortality outcomes
│   └── adjusted/                   adjusted event studies: tex/ figures/ csv/  ← analysis/parallel_trends/event_study_mortality_adjusted.R (Rishad)
└── crown_results.tex               stand-alone LaTeX doc of crown_did/adjusted + crown_event_studies/adjusted
```

Inside an event-study folder: `tex/` holds the LaTeX tables and figure blocks, `figures/` the graphs (`.pdf` for LaTeX, `.png` to look at), `csv/` every number.

## Overleaf

The Crown event-study LaTeX files refer to their figures and to each other by their path from this folder (`crown_event_studies/adjusted/...`), and their labels start with `crown_`, so they can sit in the same Overleaf project as the mortality tables. To add them, upload the `crown_event_studies/` folder (and `crown_did/` for the adjusted DiD) to the Overleaf project root, keeping the folder structure, then add:

```latex
\subsection{Crown Event Studies (economic outcomes)}
\input{crown_event_studies/adjusted/tex/es_adjusted_combined_all}
\input{crown_event_studies/adjusted/tex/figures_combined}
```

for the adjusted DiD:

```latex
\input{crown_did/adjusted/did_adjusted_combined}
\input{crown_did/adjusted/did_adjusted_all_pepfar}
\input{crown_did/adjusted/did_adjusted_cop}
```

and, for the outcome-transformation robustness check:

```latex
\input{crown_did/transformations/did_transformations_all_pepfar}
\input{crown_did/transformations/did_transformations_cop}
\input{crown_did/transformations/pretrends_transformations_all_pepfar}
\input{crown_did/transformations/pretrends_transformations_cop}
```

and, for the linear-time-trend check:

```latex
\input{crown_event_studies/robustness/tex/linear_trends_all_pepfar}
\input{crown_event_studies/robustness/tex/linear_trends_cop}
```
