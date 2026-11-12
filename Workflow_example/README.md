# Workflow example

This example uses the microscopy eukaryotic phytoplankton dataset at the genus level. It constructs a taxonomy-derived proxy tree and calculates pairwise pMST, betaNTI and RCbray. No plotting is included.

## Inputs

- `IdenPlantGenu_Asm.csv`: abundance table used in assembly analysis; 56 genera in rows and 48 mesocosms in columns.
- `IdenPlantGenu_Asm_Taxonomy.csv`: Phylum, Class, Order, Family and Genus assignments for these 56 genera.
- `treatment.csv`: sample names in the first column and treatment groups in the `Treatment` column.

The abundance and grouping tables are unchanged copies of the released inputs. The taxonomy table is an unchanged copy of the existing microscopy eukaryotic phytoplankton genus-level classification table, with five columns: Phylum, Class, Order, Family and Genus. No classification entries were reconstructed from a tree or supplemented from sequencing annotations.

## Run

Install NST, iCAMP, vegan and ape, then set this folder as the R working directory and run:

```r
source("Workflow.R", encoding = "UTF-8")
```

The script builds the proxy tree using `ape::as.phylo(..., collapse = FALSE)` and assigns a length of 1 to every branch. It removes zero-abundance samples and taxa, matches the abundance table to the tree, and converts each retained sample to relative abundance. All 48 mesocosms have nonzero total abundance, giving 120 within-treatment sample pairs across eight treatments.

For genus-level analysis, the taxonomy columns and tree formula include Phylum, Class, Order, Family and Genus. For species-level analysis, include the Species column when converting the taxonomy table to factors and extend the formula to `~Phylum/Class/Order/Family/Genus/Species`. Use the corresponding species-level abundance and taxonomy files.

The calculation settings match the released analysis scripts: 1,000 randomizations, one worker, seed 123, and abundance-weighted bMNTD with tip-label shuffling for pNST. One pNST run returns both pairwise pMST (`MST.ij.bMNTD`) and betaNTI (`bNTI.wt`). RCbray is calculated independently with `iCAMP::RC.cm`; the RC statistic returned internally by NST is not used as RCbray.

## Outputs

- `IdenPlantGenu_Asm_ProxyTree.tre`
- `IdenPlantGenu_Asm_ProxyTree_pMST_pairwise.csv`
- `IdenPlantGenu_Asm_ProxyTree_betaNTI_RCbray_pairwise.csv`

Outputs are saved in the working directory and overwritten on rerunning. No relative abundance matrices, distance files or plots are saved. RCbray is not guaranteed to be identical across repeated runs because of the package's internal random-number implementation.
