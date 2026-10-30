# Metagenomic assessment of airborne and surface viral signals and worker acceptability during a partial-lid intervention in a Cambodian live bird market

## Overview
This repository contains the analysis code for the study:

**"Metagenomic assessment of airborne and surface viral signals and worker acceptability during a partial-lid intervention in a Cambodian live bird market"**

Mechanical defeathering can disperse viral material in live poultry markets, but machine modifications must fit routine work. This nonrandomized field study evaluated partial lids on defeathering machines at Orussey Market, Cambodia, across three sampling days with lids and three without.

Targeted RT-qPCR and target-enriched metagenomics characterised viral RNA in air and on nearby surfaces, with surface samples collected before and after a three-hour sampling period. Worker questionnaires assessed usability and acceptability. The environmental findings were interpreted alongside workers' experiences to identify priorities for further evaluation.


## Key Findings

- Sequencing detected 73 poultry-associated viral taxa.
- Airborne richness was descriptively lower on days with lids, but targeted viral detections did not change consistently, and composition profiles overlapped.
- Surface profiles changed from the pre-sampling to post-sampling period on days with and without lids; the study did not establish a lid-associated difference in these changes.
- Three of five workers felt safer, but all reported discomfort and perceived slower work. None wanted to use the lid in its tested form; three would continue if it were permanently integrated into the machine.
- With only three independently allocated days per group, the environmental effect remains uncertain.

---

## Data Sources
- Air and nearby surface samples collected at Orussey Market, Cambodia, across three sampling days with partial lids and three without, including surface samples before and after a three-hour sampling period.
- Targeted RT-qPCR and target-enriched metagenomic characterisation of viral RNA; poultry-associated taxa were assessed in the environmental metagenomic analyses.
- Questionnaires assessing five workers' experiences of safety, comfort, work speed and willingness to continue using the lid.
- The script requires five separately supplied workbooks: `Original_Virus_Count.xlsx`, `Species_Count.xlsx`, `Species_RPM.xlsx`, `Metadata.xlsx` and `Avian_Species_List.xlsx`. Data files are not included with the script and README.

---

## Repository Contents
- `LBM_Intervention.R`: environmental metagenomic analysis script. It includes optional count-report preparation functions, negative-control diagnostics, alpha diversity, poultry–air Simpson diversity differences, observed-taxon plots, Bray–Curtis and quantitative Jaccard ordination, exploratory PERMANOVA, distances to the poultry median in the first two ordination axes, edgeR differential abundance, relative-abundance plots and a rank-scaled RPM heatmap. Targeted RT-qPCR and worker questionnaire analyses are not implemented in this script.


---

## Implications
The study identifies environmental and operational questions that should guide further evaluation of practical controls for poultry processing in markets. The environmental effect of the partial lid remains uncertain, and workers' reports of discomfort and slower work indicate that design and integration into routine processing require further attention.

## Citation
If you use this repository in your research, please cite:

**[Your Full Citation Here]**
