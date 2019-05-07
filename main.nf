#!/usr/bin/env nextflow

// help
def helpMSG() {
    log.info """

    Usage:

    nextflow run main.nf --input Accession_list.txt --cpus 8

    --input       a list of accession numbers, one accession number per line, no headers
                    e.g. do a 'cut -f2' on a blastn query with '-outfmt6'

    Options:
    --cpus        max cores [default $params.cpus]

    Results are stored in cluster_results/
    """.stripIndent()
}

if (params.help) { exit 0, helpMSG() }
if (params.input == '') {exit 1, "--input is a required parameter, see --help"}

Channel.fromPath(params.input)
        .splitCsv()
        .unique()
        .map { it }
        .set { accession_list }

//accession_list.subscribe { println "Got: ${it}" }

process download {
    maxForks 2
  input:
    val(line) from accession_list
  output:
    set val("${line[0]}"), file("*.gb") into download_filter_ch
  shell:
    """
    efetch -db nuccore -id !{line[0]} -format gb > !{line[0]}.gb
    while ! [ -s !{line[0]}.gb ]; do
      efetch -db nuccore -id !{line[0]} -format gb > !{line[0]}.gb
    sleep 1
    done
    """
}

process filter {
    // filtering based on circular (complete sequence) and a plasmid definition
    publishDir "${params.output}/genbank_files_filtered", mode: 'copy', pattern: "${name}_filtered.gb"
  input:
    set val(name), file(file) from download_filter_ch
  output:
    set val(name), file("${name}_filtered.gb") optional true into filter_getFasta_ch
  script:
    """
    if head -1 ${file} | grep -q "circular BCT"; then
      if grep -q -m1 "/plasmid=" ${file}; then
        cat ${file} > ${name}_filtered.gb
    fi fi
    """
}

process getFasta {
    publishDir "${params.output}/fasta_files_filtered/", mode: 'copy', pattern: "*.fasta"
    errorStrategy 'ignore'
  input:
    set val(name), file(genbank) from filter_getFasta_ch
  output:
    file("${name}.fasta") optional true into getFasta_clustering_ch
  script:
    """
    seqret -sequence ${genbank} -sformat gb -outseq ${name}.fasta
    """
}

process clustering {
    cpus = "${params.cpus}"
    publishDir "${params.output}", mode: 'copy', pattern: "output.clstr"
    publishDir "${params.output}", mode: 'copy', pattern: "representatives.txt"
    publishDir "${params.output}", mode: 'copy', pattern: "output.fasta"
  input:
    file '*.fasta' from getFasta_clustering_ch.collect()
  output:
    file("output.clstr") into clustering_createChannel
    file("representatives.txt") into name_list
    file("output.fasta") into clustering_multifastafile_ch
  script:
    """
    cat *.fasta > all.fasta
    psi-cd-hit.pl -i all.fasta -o output -aL 0.7 -aS 0.7 -circle 1 -prog megablast -para ${params.cpus}
    psi-cd-hit.pl -i all.fasta -o output -prog blastn -para ${params.cpus} -c 0.6 -s "-evalue 10E-100 -max_target_seqs 100000 -qcov_hsp_perc 10 -max_hsps 10"
    grep ">" output > representatives.txt
    cp output output.fasta
    """
}

/*
#psi-cd-hit.pl -i all.fasta -o output_blast_mods -prog blastn -para 8 -c 0.3 -s "-evalue 10E-100 -max_target_seqs 100000 -qcov_hsp_perc 10 -max_hsps 10"
- cluster 100, small cluster test bestanden
# problem sind die langen representatives im cluster
-> entweder -c0.3 langsam erhöhen
-> oder aL benutzen
#"
-> das ist noch wessentlich besser. mach damit mal eine analyse mit 1.) sourmash 2.) kompletten workflow um zu sehen wie die tabelle aussieht

*/

process save_representatives {
    publishDir "${params.output}/representatives/${name}", mode: 'copy', pattern: "*.fasta"
  input:
    set val(id), val(fasta) from clustering_multifastafile_ch.splitFasta( record: [id: true, seqString: true ])
                                                         .map { record -> tuple(record.id, record.seqString)}
  output:
    set val(id), file("*.fasta") into (representatives_prokka_ch, representatives_abricate_ch)
  script:
    """
    printf ">${id}\n${fasta}" > ${id}.fasta
    """
}

process prokka {
    publishDir "${params.output}/representatives/${name}", mode: 'copy', pattern: "prokka.gbk"
  input:
    set val(name), file(fasta) from representatives_prokka_ch
  output:
	  set val(name), file("prokka.gbk") into annotation_done_ch
  shell:
    """
    prokka --outdir output --prefix annotation !{fasta}
    mv output/annotation.gbk prokka.gbk
    """
}

process abricate {
    publishDir "${params.output}/representatives/${name}", mode: 'copy', pattern: "*.tab"
  input:
    set val(name), file(fasta) from representatives_abricate_ch
  output:
	  set val(name), file("*.tab") into annotation2_done_ch
  shell:
    """
  	abricate !{fasta} --nopath --quiet --mincov 80 --db ncbi > ncbi.tab
    abricate !{fasta} --nopath --quiet --mincov 95 --db plasmidfinder > plasmidfinder.tab
    abricate !{fasta} --nopath --quiet --mincov 35 --db transposon > transposon.tab
    """
}
