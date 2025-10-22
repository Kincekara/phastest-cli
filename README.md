# phastest-cli
[![docker size](https://badgen.net/docker/size/kincekara/phastest-cli/0.1?icon=docker&label=Docker)](https://hub.docker.com/r/kincekara/phastest-cli/tags)

This work is still under development!

Command line version of the [PHASTEST](https://phastest.ca/) web app, which is designed to support the rapid identification, and annotation of prophage sequences within bacterial genomes and plasmids.

Some of the scripts have been modified to improve portability and standalone functions.

## Installation
Docker image with lite database is available at [Dockerhub](https://hub.docker.com/r/kincekara/phastest-cli)

All databases can be found at https://phastest.ca/databases

## Quickstart
```bash
# GenBank accession
phastest -i genbank -a [accession_number]

# Contigs
phastest -i contig -s /path/to/contigs.fa

# Fasta file
phastest -i fasta -s /path/to/seq.fna
```

## Citation
Wishart, D. S., Han, S., Saha, S., Oler, E., Peters, H., Grant, J., Stothard, P., & Gautam, V. (2023). PHASTEST: faster than PHASTER, better than PHAST. Nucleic Acids Research, 51(W1), W443–W450. https://doi.org/10.1093/nar/gkad382

Wishart, D. S., Han, S., Saha, S., Oler, E., Peters, H., Grant, J., Stothard, P., & Gautam, V. (2023). PHASTEST: faster than PHASTER, better than PHAST. Nucleic Acids Research, 51(W1), W443–W450. https://doi.org/10.1093/nar/gkad382

‌Zhou, Y., Liang, Y., Lynch, K. H., Dennis, J. J., & Wishart, D. S. (2011). PHAST: A Fast Phage Search Tool. Nucleic Acids Research, 39(suppl), W347–W352. https://doi.org/10.1093/nar/gkr485
‌

## License
[![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg

