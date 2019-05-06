# KPC Plasmid analysis

![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

* Nextflow Workflow using Docker

# Dependencies

* Nextflow
* Docker

# Installation

```bash
# install nextflow
curl -s https://get.nextflow.io | bash
mv nextflow $HOME/bin
# simple install Docker:
  # alternatively see https://docs.docker.com/install/ to use the current version
apt install docker-ce
```

# Execution

```bash
git clone replikation/plasmid_analysis && cd plasmid_analysis
nextflow run main.nf --input /Accessionlist/Accessionlist_KPC_nt_archive_all_Mai_2019.txt
```
* `Accessionlist_KPC_nt_archive_all_Mai_2019.txt` is given as a example feel free to use another Accession list
  * the list has to contain one accession per line
  * no headers


* a simplified result collection can be used via `bash collect_results.sh`
