# KPC Plasmid analysis

![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

* KPC plasmid clustering workflow with Nextflow
* workflow is only tested and intended for KPC-plasmids

# Installation

## Dependencies

* [Nextflow](https://www.nextflow.io/index.html)
* [Docker](https://docs.docker.com/install/)

**A example installation for Docker and Nextflow:**

```bash
# NEXTFLOW
# java runtime for nextflow, test if installed via 'java -version'
  apt install openjdk-11-jre
# install nextflow
  curl -s https://get.nextflow.io | bash
# move 'nextflow' to a $PATH location, or add 'nextflow' to $PATH
  mv nextflow /usr/local/bin

# DOCKER
# Install via https://docs.docker.com/install/ (most recent) or:
  apt install docker-ce
# Add docker group to $USER
  sudo usermod -a -G docker $USER
# restart terminal or reboot if neccessary
```

# Usage


````bash
nextflow run replikation/plasmid_analysis --input list_of_accession_numbers.txt
````

* reproduce the KPC clustering

```bash
# clone the git and navigate into
git clone https://github.com/replikation/plasmid_analysis.git
cd plasmid_analysis
# execute Nextflow using the Accessionlist of May
nextflow run main.nf --input Accessionlist/Accessionlist_KPC_nt_archive_all_May_2019.txt --cpus 8      
```

* `Accessionlist_KPC_nt_archive_all_May_2019.txt` is given as a example
* feel free to use another accession list:
  * only one accession number per line
  * no headers
* results are stored in `cluster_results`

+ a simplified result collection for KPC-plasmids can be used via `bash collect_results.sh`

## Resume

* nextflow can resume a run
* add `-resume` to the last nextflow command

## help message

* `nextflow run main.nf --help`

```
Usage:

nextflow run main.nf --input Accession_list.txt

--input       a list of accession numbers, one accession number per line, no headers
                e.g. do a 'cut -f2' on a blastn query with '-outfmt6'

Options:
--cpus        max cores [default 8]

Results are stored in cluster_results/
```

# Additional Information
* Dockerfiles created for this workflow are located under `Dockerfiles/`
* the `work` folder can be removed afterwards

## prokka
* prokka uses a ncbi dependency `tbl2asn` which may expire at some point [see here](https://github.com/tseemann/prokka/issues/139) or [here](https://github.com/tseemann/prokka/issues/215)
  * this workflow will use prokka:latest instead to avoid this issue in the future
