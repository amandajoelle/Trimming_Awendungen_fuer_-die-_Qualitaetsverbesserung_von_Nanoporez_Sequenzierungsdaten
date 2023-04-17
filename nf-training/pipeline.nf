nextflow.enable.dsl=2

/*if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}*/

params.publishDir = './results'

params.input = "$projectDir/data/test.fastq"




singularity_prinseq = 'https://depot.galaxyproject.org/singularity/prinseq:0.20.4--hdfd78af_5'
singularity_nanofilt = 'https://depot.galaxyproject.org/singularity/nanofilt:1.1.3--py35_0'
singularity_flye = 'https://depot.galaxyproject.org/singularity/flye:2.9--py310h590eda1_1'
singularity_quast = "https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1"



// Étape de nettoyage avec prinseq
process prinseq {
  container = singularity_prinseq
  publishDir "${params.publishDir}/prinseq", mode: 'copy'
  
  input:
    path inputs
  output:
    file 'trimmed.fastq'
  script: 
    """
    prinseq-lite.pl -fastq ${inputs} -out_good trimmed -trim_left 20
    """
}

/*process convert {
  input:
    path('trimmed.fastq')
  output:
    file('trimmed.fasta')
  script:
    """
   sed -n '1~4s/^@/>/p;2~4p'  ${params.input} > trimmed.fasta
    """
}*/

// Étape de nettoyage avec nanofilt

process nanofilt {
    container = singularity_nanofilt
    publishDir "${params.publishDir}/nanofilt", mode: 'copy'

    input:
        path inputs

    output:
        file 'filtered.fastq'

    script:
        """
        NanoFilt -q 10 < ${inputs} > filtered.fastq
        """
}


// Étape d'assemblage du fichier prinseq avec flye
process flye_prinseq {

  publishDir "${params.publishDir}/flye_prinseq", mode: 'copy'
  //params.input = "$projectDir/results/prinseq/trimmed.fastq" 
  container = singularity_flye
  input:
    path (prinseq_out)
  output:
    file 'assembly_prinseq.fasta'
  script:
    """
    flye --nano-raw ${prinseq_out} -o assembly
    mv assembly/assembly.fasta assembly_prinseq.fasta
    """
}

// Étape d'assemblage du fichier nanofilt avec flye
process flye_nanofilt {
  container = singularity_flye
  input:
    path(nanofilt_out) 
  output:
    file('assembly_nanofilt.fasta')
  script:
    """
    flye --nano-raw ${nanofilt_out} -o assembly
    mv assembly/assembly.fasta assembly_nanofilt.fasta
    """
}

// Étape d'assemblage du fichier fastq en imput avec flye
process flye_origin {
  publishDir "${params.publishDir}/flye_origin", mode: 'copy'

  container = singularity_flye

  input:
    path input
    path prinseq_ch
    path nanofilt_ch

  output:
    file 'assembly_input.fasta'
    file 'assembly_prinseq_ch.fasta'
    file 'assembly_nanofilt_ch.fasta'
  script:
    """
    flye --nano-raw ${input} -o assembly_input
    flye --nano-raw ${prinseq_ch} -o assembly_prinseq_ch
    flye --nano-raw ${nanofilt_ch} -o assembly_nanofilt_ch
    mv assembly_input/assembly.fasta assembly_input.fasta
    mv assembly_prinseq_ch/assembly.fasta assembly_prinseq_ch.fasta
    mv assembly_nanofilt_ch/assembly.fasta assembly_nanofilt_ch.fasta
    """
}

// Étape d'évaluation des fichiers avec quast
process compare {
  publishDir "${params.publishDir}/compare", mode: 'copy'
  container = singularity_quast

  input:
    path x 
    path y 
    path z

  output:
    path "report_dir"

  script:
    """
    quast.py  -o report_dir ${x} ${y} ${z} 
    """
}

inputs_ch = Channel.fromPath(params.input) 

workflow {
  prinseq_ch = prinseq(inputs_ch)
  nanofilt_ch = nanofilt(inputs_ch)
  flye_origin_ch = flye_origin(inputs_ch, prinseq_ch, nanofilt_ch)
  compare(flye_origin_ch)
}