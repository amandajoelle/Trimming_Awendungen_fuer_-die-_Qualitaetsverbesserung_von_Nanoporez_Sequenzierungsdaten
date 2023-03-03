#!/usr/bin/env nextflow
if (params.imput1) {
  params.input1 = params.imput1
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput1"
  System.exit(1)
}

if (params.imput2) {
  params.input2 = params.imput2
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput2"
  System.exit(1)
}


//if (params.imput4) {
// params.input4 = params.imput4
//} else {
//println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput4"
//System.exit(1)
//}

singularityImage = "https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1"

process compare {
  container = singularityImage

  input:
    path params.input1
    path params.input2
    //path params.input3
    //path params.input4

  output:
    path "report_dir"

  script:
    """
    quast.py  -o report_dir ${params.input1} ${params.input2} 
    """
}

workflow {
  compare(params.input1, params.input2).view()
}