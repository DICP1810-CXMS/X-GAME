In cross-linking mass spectrometry (XL-MS) data analysis, each spectrum typically identifies a pair of cross-linked peptides. 
However, due to peptide sequence sharing across multiple proteins, each peptide sequence may correspond to multiple potential cross-linking sites. 
While this phenomenon also exists in conventional proteomics, the complexity is significantly amplified in XL-MS data - 
the combinatorial pairing of two peptides creates substantially more ambiguous matches between cross-linked peptides and protein interaction sites. 
This leads to the inclusion of numerous ambiguous protein-protein interactions (PPIs) and site pairs in current analytical pipelines, 
which severely compromises the reliability of subsequent protein interaction network construction and structural analysis.

To address this critical issue, we developed the "stringent_filtering_of_crosslinks" analytical tool. 
This tool implements a rigorous filtering strategy that retains only the most confident cross-linking site pair for each spectrum when reprocessing [pLink3](http://pfind.net/software/pLink/index.html) search results.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1) stringent filtering of cross-links.py
The program processes cross-linking mass spectrometry (XL-MS) data from pLink2/pLink3 search results to:
✅ Separate intra-protein (within-protein) and inter-protein (protein-protein interaction, PPI) cross-links.
✅ Reduce ambiguous matches by retaining only the highest-confidence cross-link site pair per spectrum.
✅ Annotate protein pairs with gene names using UniProt-Gene mapping.
✅ Generate structured Excel outputs for downstream analysis.

#####Input

CSV Files:
pLink3 output files (.csv) containing cross-linked peptide identifications.
Required columns:
Proteins, Peptide, Modifications, Protein_Type (Intra/Inter-Protein/None).
UniProt-Gene Mapping File:
Excel file (e.g., uniprot-gene-20220927.xlsx) with columns:
Entry (UniProt ID), Gene (Gene symbol)

#####Running the program
To run the program, use the following command:

python3 stringent filtering of cross-links.py

The program processes pLink3 search results (.csv files) by first consolidating all identified cross-linked spectra.
It then classifies each spectrum into either intra-protein or inter-protein categories based on linkage type. 
For inter-protein spectra, the program constructs a frequency database of all potential protein-protein interactions (PPIs) and implements a stringent filtering protocol: 
when multiple candidate cross-links are identified for a single spectrum, only the highest-frequency interaction site pair is retained.
This frequency-based prioritization effectively eliminates ambiguous matches while preserving the most biologically relevant interactions, ensuring high-confidence results for downstream structural and interactome analyses.

#####Output

The program generates one Excel file with 8 sheets:

intra_protein_peptides:	Unique intra-protein cross-linked peptides (no redundancy).
intra_protein_peptides_to_sites:	Matched intra-protein peptides to site pairs with identification frequency.
intra_protein_sites:	Non-redundant intra-protein site pairs (aggregated spectra counts).
intra_proteins:	Intra-protein statistics (spectra/site counts per protein).
inter_protein_peptides:	Unique inter-protein cross-linked peptides (no redundancy).
inter_protein_peptides_to_sites:	High-confidence inter-protein peptides → site pairs (retains only the top-frequency match per spectrum).
inter_protein_sites:	Non-redundant inter-protein site pairs (aggregated spectra counts).
inter_protein_PPI	PPI statistics: (spectra/site counts per protein pair).

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
2) PPI_evaluation.R - Protein-Protein Interaction Biological Relevance Assessment
This script utilizes the [GOSemSim](https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html) R package to evaluate the biological relevance of filtered protein-protein interactions (PPIs)
by calculating Gene Ontology (GO) semantic similarity across three domains:
Cellular Component (CC)
Molecular Function (MF)
Biological Process (BP)

#####Inuput
PPI result file (.csv) containing the following columns:
Protein1, Protein2, Spectra_Count, Site_Count, Gene1, Gene2

#####Running the program
Open PPI_evaluation.R in RStudio (R ≥4.4.1)
Run the script to compute GO semantic similarity scores

#####Output
Augments the input PPI table with three new columns:
CC: GO Cellular Component similarity
MF: GO Molecular Function similarity
BP: GO Biological Process similarity

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
The "data" folder contains relevant example files
