# Core network-method citation audit for rev12

## Scope

Four references were added to support the manuscript's positioning on network degree bias, normalization, benchmarking and the evidential standard for drug-repurposing validation. Bibliographic identity and claim support were checked against publisher or PubMed records on 30 July 2026.

| Manuscript ref. | Bibliographic identity | DOI | Locator supporting the manuscript claim | Verdict |
|---|---|---|---|---|
| 8 | Erten S, Bebek G, Ewing RM, Koyutürk M. *DADA: Degree-Aware Algorithms for Network-Based Disease Gene Prioritization*. BioData Mining. 2011;4:19. | `10.1186/1756-0381-4-19` | PubMed abstract and article description identify degree-aware disease-gene prioritization and the influence of PPI topology. | supports |
| 9 | Biran H, Kupiec M, Sharan R. *Comparative Analysis of Normalization Methods for Network Propagation*. Frontiers in Genetics. 2019;10:4. | `10.3389/fgene.2019.00004` | Introduction states that raw propagation is affected by seed and candidate degree; Methods/Results compare degree-based and symmetric adjacency normalization and degree-preserving reference schemes. | supports |
| 10 | Picart-Armada S, Barrett SJ, Willé DR, Perera-Lluna A, Gutteridge A, Dessailly BH. *Benchmarking Network Propagation Methods for Disease Gene Identification*. PLoS Computational Biology. 2019;15:e1007276. | `10.1371/journal.pcbi.1007276` | Abstract and benchmarking design compare 12 algorithms, two biological networks, multiple metrics and input definitions; conclusions emphasize the impact of validation strategy and seed definition. | supports |
| 11 | Cheng F, Desai RJ, Handy DE, Wang R, Schneeweiss S, Barabási A-L, Loscalzo J. *Network-Based Approach to Prediction and Population-Based Validation of In Silico Drug Repurposing*. Nature Communications. 2018;9:2691. | `10.1038/s41467-018-05116-5` | Abstract reports interactome predictions followed by patient-level pharmacoepidemiology and mechanistic in-vitro testing; used only as an example of validation beyond ranking. | supports |

## Source links

- Ref. 8: https://pubmed.ncbi.nlm.nih.gov/21699738/
- Ref. 9: https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2019.00004/full
- Ref. 10: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007276
- Ref. 11: https://www.nature.com/articles/s41467-018-05116-5

## Claim boundary

These papers support known generic risks and validation standards. They do not establish the ACC-specific seed-exposure result. The manuscript therefore states that its increment is tracing degree, normalization and seed-definition risks into a downstream drug ranking with a fixed pharmacogenomic association map, not inventing a new propagation algorithm or claiming the first observation of network degree bias.
