
params.input = "$projectDir/data/test.fastq"

/*if (params.imput) {
  params.input = params.imput
} else {
  println "Erreur : Veuillez spécifier un fichier fastq en utilisant l'option --imput"
  System.exit(1)
}*/



// Définition du conteneur Singularity à utiliser

    singularity = 'https://depot.galaxyproject.org/singularity/nanofilt:1.1.3--py35_0'


// Définition de la tâche pour exécuter NanoFilt
process run_nanofilt {
    container = singularity
    input:
        path(params.input)

    output:
        file('filtered.fastq')

    script:
        """
        NanoFilt -q 10 < ${params.input} > filtered.fastq
        """
}

// Exécution du pipeline
workflow {
    run_nanofilt(params.input).view() 
}