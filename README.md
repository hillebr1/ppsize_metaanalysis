# Meta-analysis of phytoplankton cell size 

This repository contains raw data and code for the meta-analysis of functional consequences of cell size for pyhtoplankton performance. It is part of the paper
Cell size as driver and sentinel of phytoplankton community structure and functioning
by 
Helmut Hillebrand, Esteban Acevedo-Trejos, Stefanie D. Moorthi, Alexey Ryabov, Maren Striebel, Patrick Thomas, Marie-Luise Schneider
to be published in Functional Ecology.

Please refer to the article regarding background, search and inclusion criteria for the literature, and interpretation of the results.

## Data and Scripts

The analyses are coded in "code_sizeMA_hh.R", which uses two data files, studies.csv and data_incl.csv. The former is the output of the systematic literature review, the latter includes the digitized data. 

### Metadata for "studies.csv"

Studies.csv is mostly identical to Table S1 in the Supporting Material of the paper. For each paper found by the search in Web of Science we give 

Paper: numbered in order of appearance from the search 

Authors: Author names

Year: year of publication

Citation: Journal, volume, pages and/or DOI

Included: whether the paper was included in our final database

System: M = marine F = freshwater, FW = Both, none = not specified

Type: type of study (M= model, E = experiment, O= observation, R= review or a combination thereof)

Level: Organizational level at which data were measured (SS = single species, MS = multiple single species, or C = communities)

Driver.Response: whether cell size appeared as driver (D) or response (R). 

Aim: Allocation to the aims of our review as stated in the paper.

Part.of.meta-analysis: Whether data were digitized for meta-analysis of Aim 1, if so, unique identifier for each each paper included. 

### Metadata for "data.csv"

MA_ID: unique identifier of study, identical to "studies.csv"

caseID: unique identifier for case within study, different cases e.g. comprising different response variables measured in the same study

system: M = marine F = freshwater, FW = Both

Phylum: Phylum the measured algal species belongs to

resp.categ: Allocation to the four response categories (as depicted in Fig. 2a) 

resp.scale: whether measurement was per cell (absolute) or specific per biomass or biovolume. 

resp.subset: subsets of measured responses withing categories, these are used to facet Fig. 3-5

res.unit: exact unit the response variables were measured in

log10cellsize: cell volume in µm³, log 10 transformed

newresp: response value in log 10 transformed units (with exception see below)

resptype: identifies variables (growth rate) that were not log-transformed
