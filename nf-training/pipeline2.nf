nextflow.enable.dsl=2

/*if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}*/
 
params.publishDir = './results'

params.input = "$projectDir/data/test.fastq"

params.input1 = "$projectDir/data/test1.fastq"




singularity_prinseq = 'https://depot.galaxyproject.org/singularity/prinseq:0.20.4--hdfd78af_5'
singularity_nanofilt = 'https://depot.galaxyproject.org/singularity/nanofilt:1.1.3--py35_0'
singularity_flye = 'https://depot.galaxyproject.org/singularity/flye:2.9--py310h590eda1_1'
singularity_quast = "https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1"


// Étape de nettoyage avec prinseq
process prinseq {
  container = singularity_prinseq
  publishDir "${params.publishDir}/prinseq", mode: 'copy'
  
  input:
    path input
  output:
    path "trimmed_${input.baseName}.fastq"
  script: 
    """
    prinseq-lite.pl -fastq ${input} -out_good trimmed_${input.baseName} -trim_left 20
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
        path input

    output:
        file "filtered_${input.baseName}.fastq"

    script:
        """
        NanoFilt -q 10 < ${input} > filtered_${input.baseName}.fastq
        """
}


// Étape d'assemblage du fichier fastq en imput avec flye
process flye_origin {

 

  publishDir "${params.publishDir}/flye_origin", mode: 'copy'

  container = singularity_flye

  input:
    path input

  output:
    file "assembly_${input.baseName}.fasta"
    
  script:
    """
    flye --nano-raw ${input} -o assembly_${input.baseName}
   
    mv assembly_${input.baseName}/assembly.fasta assembly_${input.baseName}.fasta
  
    """

 
}

// Étape d'évaluation des fichiers avec quast
process compare {
  publishDir "${params.publishDir}/compare", mode: 'copy'
  container = singularity_quast

  input:
    path x 

  output:
    path "report_dir"

  script:
    """
    quast.py  -o report_dir ${x}
    """
}



workflow {
  inputs_ch = Channel.fromPath(params.input) 
  prinseq_ch = prinseq(inputs_ch)
  nanofilt_ch = nanofilt(inputs_ch)
  
  concat_ch = inputs_ch.concat(prinseq_ch, nanofilt_ch)
  //inputs_ch.concat(prinseq_ch, nanofilt_ch).set {inputs_ch}
  
  //concat_ch = inputs_ch
  //concat_ch = concat_ch.concat(prinseq_ch, nanofilt_ch)
  //concat_ch.view{ "val : $it" }
  flye_origin_ch = flye_origin(concat_ch)
  flye_origin_ch.collect().view()
  compare(flye_origin_ch.collect())
}