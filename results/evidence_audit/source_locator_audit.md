# C4 source-locator audit

Judgment states: `supports`, `partial`, or `unsupported`. The judgment refers to the legacy drug-specific claim, not merely to whether the cited paper exists.

| Drug | Legacy claim | Source support | Locator judgment | Benchmark consequence |
|---|---|---|---|---|
| Mitotane | positive (1.00) | Guideline recommends mitotane monotherapy or EDP-mitotane for advanced ACC; evidence is not a randomized single-agent efficacy trial. | supports | Retain as strict positive candidate |
| Cabozantinib | positive (0.90) | 18 treated patients; 13/18 progression-free at four months; median PFS 6 months; no concomitant mitotane. | supports | Retain as strict positive candidate |
| Etoposide | positive (0.70) | FIRM-ACT randomized EDP plus mitotane versus streptozocin plus mitotane; no etoposide monotherapy contrast. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Doxorubicin | positive (0.70) | FIRM-ACT randomized EDP plus mitotane versus streptozocin plus mitotane; no doxorubicin monotherapy contrast. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Cisplatin | positive (0.70) | FIRM-ACT randomized EDP plus mitotane versus streptozocin plus mitotane; no cisplatin monotherapy contrast. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Gemcitabine | positive (0.65) | 145 patients; 132 received capecitabine; partial response 4.9% and stable disease 25.0%; drug-specific effect is not identifiable. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Palbociclib | positive (0.65) | Directly tested in SW-13 and NCI-H295R cells; effects differed by pRB context. | partial | Retain as preclinical context; exclude from clinical benchmark |
| Carboplatin | positive (0.60) | The cited regimen contains cisplatin rather than carboplatin. | unsupported | Remove drug-specific clinical label |
| Ribociclib | positive (0.60) | Directly tested in SW-13 and NCI-H295R cells; activity was cell-context dependent. | partial | Retain as preclinical context; exclude from clinical benchmark |
| Abemaciclib | external score only (0.55) | Mechanism-class rationale only; no abemaciclib experiment in the cited paper. | unsupported | Remove drug-specific clinical label |
| Lenvatinib | external score only (0.55) | Eight patients received lenvatinib plus pembrolizumab; two partial responses and one stable disease. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Actinomycin D | positive (0.55) | One infant received surgery plus actinomycin D cyclophosphamide and 5-fluorouracil. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Everolimus | external score only (0.50) | Direct everolimus exposure significantly inhibited growth in H295R and SW13 cells. | partial | Retain as preclinical context; exclude from clinical benchmark |
| Temsirolimus | external score only (0.50) | 26 patients received the combination; 11 had stable disease longer than six months. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Sirolimus | external score only (0.50) | The cited experiment evaluated everolimus rather than sirolimus. | unsupported | Remove drug-specific clinical label |
| Sunitinib | negative (0.35) | 35 evaluable patients; five stable disease and 24 progressive disease; median PFS 2.8 months; concomitant mitotane may reduce exposure. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |
| Osimertinib | negative (0.35) | The cited study evaluated gefitinib rather than osimertinib. | unsupported | Remove drug-specific clinical label |
| Afatinib | negative (0.25) | The cited study evaluated gefitinib rather than afatinib. | unsupported | Remove drug-specific clinical label |
| Erlotinib | negative (0.20) | Ten heavily pretreated patients received erlotinib plus gemcitabine; one minor response and eight progressed at first staging. | partial | Retain regimen/context evidence; exclude from strict binary benchmark |

All 19 records were checked against a primary paper or an ACC clinical guideline. Anthropic Claude performed a protocol-locked, blinded model-based classification of all 19 records; this was not a human rereview. After source-access completion for B02, strict eligibility and strict binary labels agreed in 19/19 records in the Claude–human traceability comparison.
