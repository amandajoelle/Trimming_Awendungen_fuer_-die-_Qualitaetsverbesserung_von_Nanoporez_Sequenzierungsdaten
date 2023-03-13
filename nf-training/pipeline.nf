#!/usr/bin/env nextflow

nextflow.enable.dsl=2

singularity_prinseq = 'https://depot.galaxyproject.org/singularity/prinseq:0.20.4--hdfd78af_5'
singularity_nanofilt = 'https://depot.galaxyproject.org/singularity/nanofilt:1.1.3--py35_0'
singularity_flye = 'https://depot.galaxyproject.org/singularity/flye:2.9--py310h590eda1_1'
singularity_quast = "https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1"

publishDir "$projectDir/data/results", mode: 'copy'
params.input = "$projectDir/data/test.fastq"

// Crée une Channel contenant le fichier d'entrée

output_dir = "output"

// Étape de nettoyage avec prinseq
process prinseq {
  container = singularity_prinseq

  input: 
    file trimme_fastq

  output: 
    file(trimmed_fastq) into prinseq_out

  script:
    """
    prinseq-lite.pl -fastq ${trimme_fastq} -out_good trimmed -trim_left 20
    """
}

// Étape de nettoyage avec nanofilt
process nanofilt {
    container = singularity_nanofilt
    input: 
        file 'filter.fastq' from input_sequences
    output: 
        file 'filtered.fastq' into nanofilt_out
    script:
        """
        NanoFilt -q 10 < filter.fastq > filtered.fastq
        """
}

// Étape d'assemblage du fichier prinseq avec flye
process flye_prinseq {
  
  container = singularity_flye
  input: 
    file x from prinseq_out
  output: 
    file 'assembly_prinseq.fasta' into flye_prinseq_out
  script:
    """
    flye --nano-raw ${prinseq_out} -o assembly
    """
}

// Étape d'assemblage du fichier nanofilt avec flye
process flye_nanofilt {
  container = singularity_flye
  input: 
    file x from nanofilt_out
  output: 
    file 'assembly_nanofilt.fasta' into flye_nanofilt_out
  script:
    """
    flye --nano-raw ${nanofilt_out} -o assembly
    mv assembly/assembly.fasta assembly.fasta
    """
}

// Étape d'assemblage du fichier fastq en imput avec flye
process flye_origin {
  publishDir params.publishDir 
  container = singularity_flye
  input: 
    file 'filter.fastq' from input_sequences
  output: 
    file 'assembly_origin.fasta' into flye_origin_out
  script:
    """
    flye --nano-raw ${input_sequences} -o assembly
    """
}

// Étape d'évaluation des fichiers avec quast
process compare {
  container = singularity_quast

  input:
    file x from flye_prinseq_out
    file y from flye_nanofilt_out
    file z from flye_origin_out
    //path ('assembly_origin.fasta')

  output:
    path "report_dir"

  script:
    """
    quast.py  -o report_dir ${x} ${y} ${z}  
    """
}

workflow {
  //  Input data is received through channels
  input_ch = Channel.fromPath(params.input)
  prinseq(input_ch)
  prinseq.out.view()
}
