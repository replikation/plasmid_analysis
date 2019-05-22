# KPC Plasmid analysis

![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

* KPC plasmid clustering workflow with Nextflow
* workflow is only tested and intended for KPC-plasmids
* could work similar for other beta-lactamase resistance plasmids (untested)

# Installation

## Dependencies

* [Nextflow](https://www.nextflow.io/index.html)
* [Docker](https://docs.docker.com/install/)

* a example installation for docker and Nextflow looks like this:

```bash
# NEXTFLOW
# install java runtime for nextflow
  # you can check if java is allready installed via 'java -version'
  sudo apt install openjdk-11-jre
# install nextflow
  curl -s https://get.nextflow.io | bash
# move 'nextflow' to a $PATH location, or add 'nextflow' to $PATH
  sudo mv nextflow /usr/local/bin

# DOCKER
# install via:
  apt install docker-ce
  # or via https://docs.docker.com/install/ (most recent)
# Add the docker group to $USER
  sudo usermod -a -G docker $USER
# restart terminal or reboot if neccessary
```

# Usage

```bash
git clone https://github.com/replikation/plasmid_analysis.git && cd plasmid_analysis
nextflow run main.nf --input Accessionlist/Accessionlist_KPC_nt_archive_all_May_2019.txt --cpus 8      
```

* `Accessionlist_KPC_nt_archive_all_May_2019.txt` is given as a example
* feel free to use another accession list:
  * only one accession number per line
  * no headers
* results are stored in `cluster_results`

+ a simplified result collection for KPC-plasmids can be used via `collect_results.sh`

## Resume

* nextflow supports a resume function add `-resume` to the full nextflow command

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
