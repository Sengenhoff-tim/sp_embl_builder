# sp_embl_builder

Builds UniProt flat files from [UniProt]([https://example.com](https://www.uniprot.org/api-documentation/uniprotkb),  [EBI](https://www.ebi.ac.uk/proteins/api/doc), and sample-specific variants. 
Builds entries for each provided accession and for each variant in the input list that can be mapped to UniProt.

If Ensembl fallback is enabled, sequences are retrieved from [Ensembl](https://rest.ensembl.org/)



## Status
s
Work in progress.

TODOS:

- handle flakiness of Ensembl (currently omitting Ensembl inference is recommended)
- make EBI request more efficient (currently per entry)
- resolve more EBI variants (variants that cannot be parsed are omitted)
- refactor for more clarity, improve separation of concerns
- fix tests
