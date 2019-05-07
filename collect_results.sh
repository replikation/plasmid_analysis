#!/usr/bin/env bash
CPU=4  # lscpu -p | egrep -v '^#' | wc -l
output_dir='cluster_results' # needs to be still changed in parallel
type xargs >/dev/null 2>&1 || { echo -e >&2 "${RED}xargs not found, please install Aborting.${NC}"; exit 1; }

echo "Result collector for the nextflow output"
echo " "
[ ! -d "$output_dir" ] && echo "Nextflow output dir '$output_dir' DOES NOT exists, please run nextflow first" && exit 1
[ -d "$output_dir" ] && echo "  $output_dir found, starting with result collection"


# parallel result collection
representatives=$(find $output_dir/representatives/*/*.fasta -type f -print0  | xargs -0 -n 1 basename -s '.fasta')
mkdir -p .collection
echo "$representatives" | xargs -I% -P ${CPU} \
  sh -c 'printf "%," > .collection/%.list.csv && \
        head -1 cluster_results/genbank_files/%*.gb | tr -s " " | cut -f3 -d " " | tr "\n" "," >> .collection/%.list.csv && \
        grep "/plasmid=" cluster_results/genbank_files/%*.gb | cut -f2 -d "\"" | tr "\n" "," >> .collection/%.list.csv && \
        if grep -q -m1 "/collection_date=" cluster_results/genbank_files/%*.gb; then
          grep "/collection_date=" cluster_results/genbank_files/%*.gb | cut -f2 -d "\"" | grep -Eo [0-9]{4} | tr "\n" "," >> .collection/%.list.csv
        else
          if grep -q % .years.csv; then
            grep % .years.csv | cut -f2 -d";" | tr "\n" "," >> .collection/%.list.csv
          else
            head -1 cluster_results/genbank_files/%*.gb | rev | cut -f1 -d "-" | rev | tr "\n" "," >> .collection/%.list.csv
          fi fi && \
        grep "ORGANISM  " cluster_results/genbank_files/%*.gb | cut -f5,6 -d " " | tr "\n" "," >> .collection/%.list.csv && \
        grep Inc cluster_results/representatives/%/plasmidfinder.tab | cut -f5 | sort | cut -f1 -d"(" | cut -f1 -d"_" | tr "\n" " " >> .collection/%.list.csv && \
          printf "," >> .collection/%.list.csv && \
        grep blaKPC cluster_results/representatives/%/ncbi.tab | cut -f5 | sort | tr -d "bla"| tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep bla cluster_results/representatives/%/ncbi.tab | cut -f5 | grep blaKPC -v | sort | tr -d "bla"| tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep bla cluster_results/representatives/%/ncbi.tab -v | tail -n+2 | cut -f5 | sort | uniq | tr " " "; " | tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        tail -n+2 cluster_results/representatives/%/transposon.tab | cut -f6,7 | tr " " "," | tr "\n" "," >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv \
        '

echo "plasmid group, Accession, length [bp], name, collection date, organism, Inc(s), KPC type(s), other beta-lactamases, other AB resistances, \
      Tn4401 coverage (1st hit), visual (1st hit), Tn4401 coverage (2nd hit), visual (2nd hit)" > plasmid_metadata.csv
i="0"
for tmpfile in .collection/*.list.csv; do
  i=$((i+1))
  num=$(printf "%03d" $i)
  printf "$num," >> plasmid_metadata.csv
  cat $tmpfile >> plasmid_metadata.csv
  echo "" >> plasmid_metadata.csv
done

echo "Results saved to 'plasmid_metadata.csv'"
rm -fr .collection/




# get all inc types
# for x in */plas* ; do echo "________________" ; tail -n+2 $x | cut -f5 | sort | grep "Inc" | cut -f1 -d"(" | cut -f1 -d"_" | uniq -c ; done
