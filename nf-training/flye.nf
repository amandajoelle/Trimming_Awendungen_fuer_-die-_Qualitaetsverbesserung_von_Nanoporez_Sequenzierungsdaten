#!/usr/bin/env nextflow

if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}


singularity = 'https://depot.galaxyproject.org/singularity/flye:2.9--py310h590eda1_1'

process assemble {
  container = singularity
  input:
    path (params.input)
  output:
    file('assembly.fasta')
  script:
    """
    flye --nano-raw ${params.input} -o assembly
    mv assembly/assembly.fasta assembly.fasta
    """
}

workflow {
  assemble(params.input).view()
}