#!/usr/bin/env nextflow



if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}

singularity = 'https://depot.galaxyproject.org/singularity/prinseq:0.20.4--hdfd78af_5'

process trim {
  container = singularity
  input:
    path(params.input)
  output:
    file('trimmed.fastq')
  script:
    """
    prinseq-lite.pl -fastq ${params.input} -out_good trimmed -trim_left 20
    """
}

process convert {
  input:
    path('trimmed.fastq')
  output:
    file('trimmed.fasta')
  script:
    """
   sed -n '1~4s/^@/>/p;2~4p'  ${params.input} > trimmed.fasta
    """
}

workflow {
  //scatter(params.input)
  trim(params.input)
  convert(trim.out).view()
}