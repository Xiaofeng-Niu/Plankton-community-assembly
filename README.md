# Plankton community assembly analysis

This project provides analysis data and R scripts for reproducing the core results of the study. Raw sequencing files are not included.

## Data

Abundance tables contain taxa in rows and samples in columns, with unique taxon names in the first column. In `treatment.csv`, the first column contains sample names and the `Treatment` column identifies treatment groups.

| Filename component | Meaning |
|---|---|
| `Iden` / `Seq` | Microscopy / amplicon sequencing |
| `Plant` / `Zoo` | Phytoplankton / zooplankton |
| `Genu` / `Spec` / `Otu` | Genus / species-level analysis table / OTU |
| `Pre` | After freshwater target-taxon selection, before taxon exclusions for assembly analysis |
| `Asm` | Eukaryotic phytoplankton or zooplankton data used in assembly analysis |
| `AsmCy` | Microscopy phytoplankton data including cyanobacteria, used in assembly analysis |
| `ProxyTree` / `SeqTree` | Taxonomy-derived proxy tree / sequence-based tree |

Filenames combine the dataset identifier, analysis stage and, for tree files, tree type, separated by underscores. Examples:

| Example filename | Description |
|---|---|
| `IdenPlantSpec_Pre.csv` | Microscopy eukaryotic phytoplankton species-level abundance table before assembly filtering |
| `IdenPlantSpec_Asm.csv` | Corresponding abundance table used in assembly analysis |
| `IdenPlantSpec_Asm_ProxyTree.tre` | Proxy tree corresponding to this abundance table |
| `IdenPlantSpec_AsmCy.csv` | Microscopy phytoplankton species-level abundance table including cyanobacteria, used in assembly analysis |
| `IdenPlantSpec_AsmCy_ProxyTree.tre` | Proxy tree corresponding to the abundance table including cyanobacteria |
| `SeqZooGenu_Asm.csv` | Amplicon sequencing zooplankton genus-level abundance table used in assembly analysis |
| `SeqZooGenu_Asm_SeqTree.nwk` | Sequence-based tree corresponding to this sequencing abundance table |
| `SeqZooGenu_Asm_ProxyTree.tre` | Proxy tree corresponding to the same sequencing abundance table |

For each dataset, `Asm` contains the taxa retained from `Pre` for assembly analysis. For microscopy phytoplankton, `AsmCy` adds the cyanobacterial taxa included in the analysis to the corresponding `Asm` taxon set. `Pre` does not denote unprocessed sequencing data. A species-level analysis table does not imply that every entry has been identified to a named species.

Match abundance tables and trees by their dataset and stage prefixes. The same sequencing abundance table can be analyzed using either a sequence-based tree or a proxy tree. The initial `Seq` identifies the data source; the final `SeqTree` or `ProxyTree` identifies the tree type. The calculation scripts use `Asm` or `AsmCy` inputs; `Pre` tables should not replace these inputs.

Output filenames append the metric to the input prefix. For example:

- `IdenPlantSpec_Asm_tMST_pairwise.csv`: pairwise tMST results for this abundance table.
- `SeqZooGenu_Asm_SeqTree_pMST_pairwise.csv`: pairwise pMST results for the sequencing abundance table using its sequence-based tree.
- `SeqZooGenu_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv`: pairwise βNTI and RCbray results for the same sequencing abundance table using its proxy tree.
- `IdenPlantSpec_Asm_ProxyTree_randomization_pairwise.csv`: pairwise results for the reference proxy tree and 999 randomized proxy trees; `replicate = 0` denotes the reference tree and replicates 1–999 denote randomized trees.

## Required R packages

The required packages and versions used for this analysis are listed below to help reproduce the software environment. Other versions may also work.

| R package | Version used | Purpose |
|---|---|---|
| NST | 3.1.10 | tMST, pMST and βNTI calculations |
| iCAMP | 1.9.1 | RCbray calculation |
| vegan | 2.6-6.1 | Relative abundance transformation |
| ape | 5.8 | Tree manipulation, taxon distances and tree randomization |
| mgcv | 1.9-1 | Classification agreement models |
| ggplot2 | 4.0.2 | Plotting |
| dplyr | 1.1.4 | Result summaries |
| tidyr | 1.3.2 | Data reshaping |
| grid | Included with R | Graphical layout and physical dimensions |
| gtable | 0.3.6 | Panel assembly |
| ragg | 1.3.2 | PNG export |
| svglite | 2.2.2 | Editable SVG export |

Install the required packages and their dependencies before running the scripts. The analysis was run in R 4.4.1. Figures use the Arial font.

## Running the scripts

Place the scripts, required abundance tables, tree files and grouping table in the same local working directory, retaining their filenames. Set this folder as the working directory in R or RStudio, then run the calculation scripts:

```r
source("tNST.R", encoding = "UTF-8")
source("pNST.R", encoding = "UTF-8")
source("bNTI_RCbray.R", encoding = "UTF-8")
source("ProxyTree_randomization.R", encoding = "UTF-8")
source("Assembly_agreement_model.R", encoding = "UTF-8")
source("Alpha_diversity.R", encoding = "UTF-8")
```

Run `Assembly_agreement_model.R` after `pNST.R`. `Alpha_diversity.R` requires only the abundance tables and can run independently. The remaining assembly calculation scripts generate their own results without requiring outputs from the other methods. Input files and calculation parameters are listed at the beginning of each script.

After the required calculations have finished, run the plotting scripts:

```r
source("Fig2.R", encoding = "UTF-8")
source("Fig3.R", encoding = "UTF-8")
source("Fig4b.R", encoding = "UTF-8")
source("Fig4d.R", encoding = "UTF-8")
source("Fig5.R", encoding = "UTF-8")
source("Fig6defg.R", encoding = "UTF-8")
source("Fig7ab.R", encoding = "UTF-8")
source("FigS2.R", encoding = "UTF-8")
source("FigS3.R", encoding = "UTF-8")
source("FigS4a.R", encoding = "UTF-8")
source("FigS5.R", encoding = "UTF-8")
```

Plotting scripts read calculated results without rerunning null models or statistical tests. To generate a particular figure, complete the calculations listed in the table below and then run its plotting script.

All outputs are saved in the working directory. Rerunning a script overwrites outputs with the same filenames. Calculations produce CSV files; plotting scripts produce SVG and PNG files with matching basenames.

## Script–figure correspondence

| Calculation script | Analysis |
|---|---|
| `tNST.R` | Pairwise tMST |
| `pNST.R` | Pairwise pMST |
| `bNTI_RCbray.R` | Pairwise βNTI and RCbray |
| `ProxyTree_randomization.R` | Sensitivity analysis using 999 randomized proxy trees |
| `Assembly_agreement_model.R` | pMST classification agreement models; requires outputs from `pNST.R` |
| `Alpha_diversity.R` | Alpha diversity indices and paired tests |

| Plotting script | Figure | Required calculations |
|---|---|---|
| `Fig2.R` | Fig. 2 | tNST, pNST, βNTI–RCbray |
| `Fig3.R` | Fig. 3 | pNST, βNTI–RCbray |
| `Fig4b.R`, `FigS4a.R` | Fig. 4b, Fig. S4a | Alpha_diversity |
| `Fig4d.R` | Fig. 4d | Assembly_agreement_model |
| `Fig5.R`, `Fig6defg.R` | Fig. 5, Fig. 6d–g | ProxyTree_randomization |
| `Fig7ab.R` | Fig. 7a, b | βNTI–RCbray |
| `FigS2.R` | Fig. S2 | tNST, pNST |
| `FigS3.R` | Fig. S3 | βNTI–RCbray |
| `FigS5.R` | Fig. S5 | tNST, pNST, βNTI–RCbray |

The model outputs `Agreement_all_model_coefficients.csv`, `Agreement_all_model_fit_summary.csv` and `Agreement_residual_intervals.csv` correspond to Tables S6, S7 and S8, respectively. Other `Agreement_*.csv` files contain model records, interaction tests, simple slopes, VIFs and record-level residuals.

## Runtime notes

The default is 1,000 internal null-model randomizations. The proxy-tree sensitivity analysis uses 999 randomized trees for each of six datasets. The full proxy-tree analysis took approximately eight hours on the computer used for this study; other calculations took several minutes to tens of minutes.

Repeated RCbray runs are not guaranteed to produce identical values because of the package's internal random-number implementation. Quantile ranges across randomized trees describe sensitivity to randomization and are not confidence intervals based on biological replicates.
