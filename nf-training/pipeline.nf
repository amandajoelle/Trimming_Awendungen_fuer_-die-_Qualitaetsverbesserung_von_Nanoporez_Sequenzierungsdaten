

//params.imput = ''
if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}


singularity_prinseq = 'https://depot.galaxyproject.org/singularity/prinseq:0.20.4--hdfd78af_5'
singularity_nanofilt = 'https://depot.galaxyproject.org/singularity/nanofilt:1.1.3--py35_0'
singularity_flye = 'https://depot.galaxyproject.org/singularity/flye:2.9--py310h590eda1_1'
singularity_quast = "https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1"

params.publishDir = './results'
params.output_prefix = 'output'

// Étape de nettoyage avec prinseq
process prinseq {
  container = singularity_prinseq
  publishDir(
        path: "${params.publishDir}/prinseq",
        mode: 'copy'
  )

  input:
    path(params.input)
  output:
    path("${params.output_prefix}_trimmed.fastq")
  
  script:
    """
    prinseq-lite.pl -fastq ${params.input} -out_good ${params.output_prefix}_trimmed -trim_left 20
    """
}


// Étape de nettoyage avec nanofilt

process nanofilt {
    container = singularity_nanofilt
    input:
        path(params.input)

    output:
        path("${params.publishDir}/${params.output_prefix}_filtered.fastq")
    
    script:
        """
        NanoFilt -q 10 < ${params.input} > ${params.output_prefix}_filtered.fastq
        """
}


// Étape d'assemblage du fichier prinseq avec flye
process flye_prinseq {
  
  container = singularity_flye
  input:
    path ("${params.publishDir}/${params.output_prefix}_trimmed.fastq")
  output:
     path ("${params.publishDir}/${params.output_prefix}_assembly_prinseq.fasta")
  
  script:
    """
    flye --nano-raw ${input} --out-dir ${params.publishDir}/output_assembly_prinseq
    """
}

// Étape d'assemblage du fichier nanofilt avec flye
process flye_nanofilt {
  container = singularity_flye
  input:
    path ("${params.publishDir}/${params.output_prefix}_filtered.fastq")
  output:
    path ("${params.publishDir}/${params.output_prefix}_assembly_nanofilt.fasta")
  publishDir params.publishDir 
  script:
    """
     flye --nano-raw ${input} --out-dir ${params.publishDir}/output_assembly_nanofilt
    """
}

// Étape d'assemblage du fichier fastq en imput avec flye
process flye_origin {
  publishDir params.publishDir 
  container = singularity_flye
  input:
    path (params.input)
  output:
    path("${params.publishDir}/${params.output_prefix}_assembly_origin.fasta")
  script:
    """
    flye --nano-raw ${input} --out-dir ${params.publishDir}/output_assembly_origin
    """
}

// Étape d'évaluation des fichiers avec quast
process compare {
  container = singularity_quast

  input:
    path ("${params.publishDir}/${params.output_prefix}_assembly_nanofilt.fasta")
    path ("${params.publishDir}/${params.output_prefix}_assembly_origin.fasta")
    path ("${params.output_prefix}/${params.output_prefix}_assembly_prinseq.fasta")
    //path ('assembly_origin.fasta')

  output:
    path "report_dir"

  script:
    """
    quast.py  -o report_dir ${input} 
    """
}


workflow {
  //scatter(params.input)
  prinseq (params.input).view()
  //convert(prinseq.out).view()
  nanofilt(params.input).view()
  flye_prinseq("${params.output_prefix}/${params.output_prefix}_trimmed.fastq").view()
  flye_nanofilt("${params.output_prefix}/${params.output_prefix}_filtered.fastq").view()
  flye_origin(params.input).view()
  compare("${params.output_prefix}/${params.output_prefix}_assembly_nanofilt.fasta", "${params.output_prefix}/${params.output_prefix}_assembly_origin.fasta", "${params.output_prefix}/${params.output_prefix}_assembly_prinseq.fasta").view()
}