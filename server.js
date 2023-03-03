const express = require('express');
const bodyParser = require('body-parser');
const { exec } = require('child_process');

const app = express();
app.use(bodyParser.urlencoded({ extended: true }));

app.post('/prinseq', (req, res) => {

    const input = req.body.input;

    const command = "nextflow run nf-training/prinseq.nf --input "+input;
    res.send(command);

    //const command = "nextflow run nf-training/prinseq.nf";

   
});

app.get('/nanofilt', (req, res) => {

    //const file = req.body.file;

    //const input = req.body.input;

    //const command = "nextflow run" +file+" --input "+input;

    const command = "nextflow run nf-training/prinseq.nf";

    exec(command, (error, stdout, stderr) => {
        if (error) {
            res.status(500).send('Une erreur s\'est produite lors de l\'exécution de la commande');
        } else {
            res.send(stdout);
        }
    });
});




app.listen(3000, () => {
    console.log('Serveur démarré sur le port 3000');
});
