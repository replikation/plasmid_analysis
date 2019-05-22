#!/usr/bin/env bash
type xargs >/dev/null 2>&1 || { echo -e >&2 "${RED}xargs not found, please install. Aborting.${NC}"; exit 1; }
#CPUs
CPU=$(lscpu -p | egrep -v '^#' | wc -l)

# Result folder
  output_dir='cluster_results' # from nextflow needs to be still changed in parallel section
  RESULT_DIR="summary"
  # Result file location
  CSV_loc="${RESULT_DIR}/plasmid_metadata.csv"
  R_loc="${RESULT_DIR}/R_sum_genetable.csv"

echo "Result collector for the nextflow output"
echo " "
[ ! -d "$output_dir" ] && echo "Nextflow output dir '$output_dir' DOES NOT exists, please run nextflow first" && exit 1
[ -d "$output_dir" ] && echo "  $output_dir found, starting with result collection"

# parallel result collection
# creates a csv table
representatives=$(find $output_dir/representatives/*/*.fasta -type f -print0  | xargs -0 -n 1 basename -s '.fasta')
mkdir -p .collection
echo "$representatives" | xargs -I% -P ${CPU} \
  sh -c 'printf "%," > .collection/%.list.csv && \
        head -1 cluster_results/genbank_files_filtered/%*.gb | tr -s " " | cut -f3 -d " " | tr "\n" "," >> .collection/%.list.csv && \
        grep "/plasmid=" cluster_results/genbank_files_filtered/%*.gb | cut -f2 -d "\"" | tr "\n" "," >> .collection/%.list.csv && \
        if grep -q -m1 "/collection_date=" cluster_results/genbank_files_filtered/%*.gb; then
          grep "/collection_date=" cluster_results/genbank_files_filtered/%*.gb | cut -f2 -d "\"" | grep -Eo [0-9]{4} | tr "\n" "," >> .collection/%.list.csv
        else
          if grep -q % .years.csv; then
            grep % .years.csv | cut -f2 -d";" | tr "\n" "," >> .collection/%.list.csv
          else
            head -1 cluster_results/genbank_files_filtered/%*.gb | rev | cut -f1 -d "-" | rev | tr "\n" "," >> .collection/%.list.csv
          fi fi && \
        grep "ORGANISM  " cluster_results/genbank_files_filtered/%*.gb | cut -f5,6 -d " " | tr "\n" "," >> .collection/%.list.csv && \
        grep Inc cluster_results/representatives/%/plasmidfinder.tab | cut -f5 | sort | cut -f1 -d"(" | cut -f1 -d"_" | tr "\n" " " >> .collection/%.list.csv && \
          printf "," >> .collection/%.list.csv && \
        grep blaKPC cluster_results/representatives/%/ncbi.tab | cut -f5 | sort | tr -d "bla"| tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep bla cluster_results/representatives/%/ncbi.tab | cut -f5 | grep blaKPC -v | sort | tr -d "bla"| tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep bla cluster_results/representatives/%/ncbi.tab -v | tail -n+2 | cut -f5 | sort | uniq | tr " " "; " | tr "\n" " " >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        grep "/country=" cluster_results/genbank_files_filtered/%*.gb | cut -f2 -d "\"" | cut -f1 -d":"| tr -d "\n" >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv && \
        tail -n+2 cluster_results/representatives/%/transposon.tab | cut -f5,6,7 | tr "\t" "," | tr "\n" "," >> .collection/%.list.csv  && \
          printf "," >> .collection/%.list.csv'

          mkdir -p $RESULT_DIR
  echo "representative, Accession, length [bp], name, collection date, organism, Inc(s), KPC type(s), other beta-lactamases, other AB resistances, country, \
        Tn type (1st hit), Tn coverage (1st hit), visual (1st hit), Tn type (2nd hit), Tn coverage (2nd hit), visual (2nd hit)" > summary/plasmid_metadata.csv
        i="0"
      for tmpfile in .collection/*.list.csv; do
        i=$((i+1))
        num=$(printf "%03d" $i)
        printf "$num," >> $CSV_loc
        cat $tmpfile >> $CSV_loc
        echo "" >> $CSV_loc
      done

echo "Results saved to '$CSV_loc'"
rm -fr .collection/

# create R overview
#create variable lists to analyse
## Incs
all_Incs=$(<$CSV_loc cut -f7 -d "," | tail -n+2 | tr " " "\n" | sort | uniq | awk 'NF')
## blas
all_blas=$(<$CSV_loc cut -f9 -d "," | tail -n+2 | tr " " "\n" | cut -f1 -d"-" | sort | uniq | awk 'NF')
## blas
all_KPCs=$(<$CSV_loc cut -f8 -d "," | tail -n+2 | tr " " "\n" | sort | uniq | awk 'NF')
## organism
all_orgs=$(<$CSV_loc cut -f6 -d "," | tail -n+2 | cut -f1 -d" " | sort | uniq | awk 'NF')

printf "gene;sample;mutated;label;year\n" > $R_loc

tail -n+2 $CSV_loc | while IFS=$',' read group Acc len name year org Inc KPC bLa oth_AB country Tn_1 Tn1_range Tn_1_cov Tn_2 Tn2_range Tn_2_cov;
    do
      export groupP=$group && export R_locP=$R_loc
      export IncP=$Inc && export bLaP=$bLa && export orgP=$org
      export kpcP=$KPC
      tranP=$(if grep -q Tn4401 <<<$Tn_1; then echo "Tn4401 or Tn4401-like"; else if grep -q Tn1721A <<<$Tn_1; then echo "Tn1721A or Tn1721A-like"; else if grep -q Tn6367 <<<$Tn_1; then echo "Tn*"; else echo "unknown" ;fi;fi;fi)
      export tranP
    # parallel inc groups
      echo "$all_Incs" | xargs -I% -P ${CPU} \
      sh -c 'if echo $IncP | grep -wq %; then printf "%;$groupP;3;Inc-group;$tranP\n" >> $R_locP; else printf "%;$groupP;0;Inc-group;$tranP\n" >> $R_locP; fi'
    # parallel beta-lactamase
      echo "$all_blas" | xargs -I% -P ${CPU} \
      sh -c 'if echo $bLaP | grep -wq %; then printf "%;$groupP;1;beta-lactamase;$tranP\n" >> $R_locP; else printf "%;$groupP;0;beta-lactamase;$tranP\n" >> $R_locP; fi'
    # parallel KPC-type
      echo "$all_KPCs" | xargs -I% -P ${CPU} \
      sh -c 'if echo $kpcP | grep -wq %; then printf "%;$groupP;1;KPC-type;$tranP\n" >> $R_locP; else printf "%;$groupP;0;KPC-type;$tranP\n" >> $R_locP; fi'
    # parallel organisms
      echo "$all_orgs" | xargs -I% -P ${CPU} \
      sh -c 'if echo $orgP | grep -wq %; then printf "%;$groupP;4;organism;$tranP\n" >> $R_locP; else printf "%;$groupP;0;organism;$tranP\n" >> $R_locP; fi'
    # features
      toxsearch=$(cat cluster_results/representatives/${Acc}/prokka.tsv)
      prokres=$(cut -f4 <<<$toxsearch)
      #umuD/C error prone DNA poly
      if grep -q "umu" <<< $prokres ; then umu=2; else umu=0; fi
      printf "umuC/D;$groupP;$umu;features;$tranP\n" >> $R_locP
      #toxin antitoxing
      if grep -q "Antitoxin" <<< $toxsearch ; then tox=2; else tox=0; fi
      printf "Toxin-antitoxin;$groupP;$tox;features;$tranP\n" >> $R_locP
      #daugther cell distro parA/B/M and stbA/B/C
      if grep -q "par" <<< $prokres ; then par=2; else par=0; fi
      if [[ "$par" == "0" ]]; then if grep -q "stb" <<< $prokres ; then par=2; fi ; fi
      printf "Plasmid-segregation;$groupP;$par;features;$tranP\n" >> $R_locP
      #SOS inhib prot psi
      if grep -q "psi" <<< $prokres ; then psi=2; else psi=0; fi
      printf "psiA/B;$groupP;$psi;features;$tranP\n" >> $R_locP
      #Mercury resistance
      if grep -w "mer[[:upper:]].." <<< $toxsearch ; then Mercuri=2; else Mercuri=0; fi
      printf "Mercuric Res.;$groupP;$Mercuri;features;$tranP\n" >> $R_locP
      #Copper resistance
      if grep -q "cop[[:upper:]].." <<< $toxsearch ; then Copper=2; else Copper=0; fi
      printf "Copper Res.;$groupP;$Copper;features;$tranP\n" >> $R_locP
      #Arsenic resistance
      if grep -q "ars[[:upper:]].." <<< $toxsearch ; then Arsenic=2; else Arsenic=0; fi
      printf "Arsenic Res.;$groupP;$Arsenic;features;$tranP\n" >> $R_locP
done
# rename CTX to CTX-M
sed -i 's/CTX/CTX-M/g' $R_locP
