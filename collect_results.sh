#!/usr/bin/env bash
CPU=4  # lscpu -p | egrep -v '^#' | wc -l
output_dir='cluster_results'
type xargs >/dev/null 2>&1 || { echo -e >&2 "${RED}xargs not found. Aborting.${NC}"; exit 1; }

echo "Result collector for the nextflow output"
echo " "
[ ! -d "$output_dir" ] && echo "Nextflow output dir '$output_dir' DOES NOT exists, please run nextflow first" && exit 1
[ -d "$output_dir" ] && echo "  $output_dir found, starting with result collection"

representatives=$(find $output_dir/representatives/*.fasta -type f -print0  | xargs -0 -n 1 basename -s '.fasta')

# parallel
mkdir -p .collection
echo "$representatives" | xargs -I% -P ${CPU} \
  sh -c 'printf "%," > .collection/%.list.csv && \
        head -1 cluster_results/genbank_files/%*.gb | tr -s " " | cut -f3 -d " " | tr "\n" "," >> .collection/%.list.csv && \
        grep "/plasmid=" cluster_results/genbank_files/%*.gb | cut -f2 -d "\"" | tr "\n" "," >> .collection/%.list.csv && \
        grep "ORGANISM  " cluster_results/genbank_files/%*.gb | cut -f5,6 -d " " | tr "\n" "," >> .collection/%.list.csv && \
        grep blaKPC cluster_results/representatives/%/ncbi.tab | cut -f5 | sort | tr -d "bla"| tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep bla cluster_results/representatives/%/ncbi.tab | cut -f5 | grep blaKPC -v | sort | tr -d "bla"| tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep Inc cluster_results/representatives/%/plasmidfinder.tab | cut -f5 | sort | cut -f1 -d"(" | cut -f1 -d"_" | tr "\n" " " >> .collection/%.list.csv && \
          printf "," >> .collection/%.list.csv && \
        grep bla cluster_results/representatives/%/ncbi.tab -v | tail -n+2 | cut -f5 | sort | tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv
        '

echo "Accession, length [bp], name, organism, KPC type(s), other beta-lactamases, other AB resistances" > plasmid_metadata.csv
for tmpfile in .collection/*.list.csv; do
  cat $tmpfile >> plasmid_metadata.csv
  echo "" >> plasmid_metadata.csv
done

echo "Results saved to 'plasmid_metadata.csv'"
rm -fr .collection/




# get all inc types
# for x in */plas* ; do echo "________________" ; tail -n+2 $x | cut -f5 | sort | grep "Inc" | cut -f1 -d"(" | cut -f1 -d"_" | uniq -c ; done
