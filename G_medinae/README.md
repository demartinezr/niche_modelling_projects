# Ecological niche and habitat availability of *Genieridium medinae* (Coleotpera: Scarabaeidae)

## Overview

This repository contains the workflow used to estimate the potential climatic distribution and available habitat of *Genieridium medinae*, an endemic dung beetle species from the Colombian Andes. The analysis was conducted to evaluate the discrepancy between the species’ potential climatic niche and the actual availability of suitable habitat, with implications for conservation assessment.

*Genieridium medinae* (Coleoptera: Scarabaeidae, Scarabaeinae) inhabits montane forests of the Central and Eastern Andes of Colombia at elevations between approximately 1800–2130 m. Dung beetles are widely recognized as indicators of ecosystem integrity and biodiversity, making this species particularly relevant for assessing the conservation status of Andean landscapes.

The purpose of this analysis is to support conservation assessment by distinguishing areas that are climatically suitable from those that currently provide appropriate habitat conditions.

## Conceptual Framework

Species distributions are shaped by multiple constraints operating at different scales. This workflow follows a hierarchical framework derived from ecological niche theory:

### Accessible area (M)

Model calibration is restricted to the area considered historically accessible to the species, based on occurrence records and dispersal assumptions. Defining M reduces model bias and prevents unrealistic extrapolation into regions that the species has likely never colonized.

### Potential climatic niche

The potential climatic niche represents geographic areas where environmental conditions are suitable for the species according to climate variables alone. This approximates the climatic component of the fundamental niche and is estimated using presence-only modeling.

### Realized habitat

Actual species distributions depend not only on climate but also on habitat structure and landscape conditions. In this workflow, habitat availability is approximated using forest cover, assuming the species is associated with intact montane forest ecosystems.

A threshold of dense canopy cover is used to represent structurally suitable habitat.

### Niche–Habitat Relationship

By intersecting climatic suitability with habitat availability, the analysis distinguishes between:

-Areas that are environmentally suitable but lack adequate habitat

-Areas that potentially support viable populations

-Regions where further field surveys may be warranted

### Conservation Context

Genieridium medinae is an endemic Andean dung beetle inhabiting mid-elevation montane forests. Members of the Scarabaeinae are widely recognized as indicators of ecosystem integrity, making them valuable for assessing landscape condition and biodiversity conservation.

Understanding the spatial relationship between climatic suitability and habitat availability is particularly important in mountain regions experiencing rapid land-use change and fragmentation.

## Intended Use

This workflow is designed for reproducible analyses in conservation biogeography and species distribution modeling, especially for:

- Habitat-dependent species

- Endemic taxa with limited distributions

- Mountain ecosystems

- Conservation planning and risk assessment

- Identification of priority survey areas

## Notes

This repository provides the modeling workflow only. Interpretation of outputs and conservation conclusions depend on additional ecological knowledge, field validation, and formal assessment criteria.
