nextflow.enable.dsl=2

/*if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}*/

params.publishDir = './results'

params.input = "$projectDir/data/test.fastq"

inputs_ch = Channel.fromPath(params.input)


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
    path (params.input)
  output:
    file 'filtered.fastq' from prinseq 
  script:
    """
    flye --nano-raw ${params.input} -o assembly
    mv assembly/assembly.fasta assembly_prinseq.fasta
    """
}

// Étape d'assemblage du fichier nanofilt avec flye
process flye_nanofilt {
  container = singularity_flye
  input:
    path ('filtered.fastq')
  output:
    file('assembly_nanofilt.fasta')
  script:
    """
    flye --nano-raw ${params.input} -o assembly
    mv assembly/assembly.fasta assembly_nanofilt.fasta
    """
}

// Étape d'assemblage du fichier fastq en imput avec flye
process flye_origin {
  publishDir "${params.publishDir}/flye_origin", mode: 'copy'

  container = singularity_flye

  input:
    path inputs
  output:
    file 'assembly_origin.fasta'
  script:
    """
    flye --nano-raw ${inputs} -o assembly
    mv assembly/assembly.fasta assembly_origin.fasta
    """
}

// Étape d'évaluation des fichiers avec quast
process compare {
  container = singularity_quast

  input:
    path ('assembly_prinseq.fasta')
    path ('assembly_nanofilt.fasta')
    //path ('assembly_origin.fasta')

  output:
    path "report_dir"

  script:
    """
    quast.py  -o report_dir ${assembly_prinseq.fasta} ${assembly_nanofilt.fasta}  
    """
}


workflow {
  //scatter(params.input)
  prinseq (inputs_ch).view()
  //convert(prinseq.out).view()
  nanofilt(inputs_ch).view()
  //flye_origin(inputs_ch).view()
  flye_prinseq(filtered.fastq).view()
  //flye_nanofilt(nanofilt.out).view()
  //compare(flye_prinseq.out, flye_nanofilt.out).view()
}