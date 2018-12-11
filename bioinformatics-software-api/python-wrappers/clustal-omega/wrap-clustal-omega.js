/**
 * Created by Vijini on 10/20/17.
 */

import cmd from 'node-cmd';
import fs from 'fs';
import uuidv4 from 'uuid/v4';
import os from 'os';
import rscript from 'js-call-r';

export default class WrapClustalOmega {
    input;

    // TODO these are the parameters
    constructor(inp = null) {
        this.input = inp;
    }

    exec() {
        
        return new Promise((resolve, reject) => {
            rscript.call(`${__dirname}/test2.R`, {
                a: this.input,
                b: 3
            })
            .then((result) => {
                console.log(result.result);
                resolve({step: 'clustal-omega', output: result.result});
            })
            .catch(err => {
                console.log('error... = ', err);
                reject("Clustal Error\n" + err);
            });
        });
/*
        // Get input file name
        const inputFilename = uuidv4();

        // Create input file with input data
        fs.writeFileSync(`${__dirname}/${inputFilename}.fasta`, this.input);

        // Get output file name
        const outputFilename = uuidv4();

        return new Promise((resolve, reject) => {
            const binary = os.type() === 'Linux' ? 'clustalo-1.2.4-Ubuntu-x86_64' : 'clustalo';
            const command = `${__dirname}/${binary} -i ${__dirname}/${inputFilename}.fasta -o ${__dirname}/${outputFilename}.fasta --auto --force -v`;

            // Execute clustal omega command
            cmd.get(command, (err, data, stderr) => {

                    if (err || stderr) {
                        reject("Clustal Error\n" + err);
                    } else {
                        // Read .fasta result
                        const content = fs.readFileSync(`${__dirname}/${outputFilename}.fasta`, 'utf8');
                        // console.log('CONTENT: ', content);
                        //console.log(content);

                        resolve({step: 'clustal-omega', output: content});
                    }

                    // Remove input file
                    fs.existsSync(`${__dirname}/${inputFilename}.fasta`) && fs.unlinkSync(`${__dirname}/${inputFilename}.fasta`);

                    // Remove output file
                    fs.existsSync(`${__dirname}/${outputFilename}.fasta`) && fs.unlinkSync(`${__dirname}/${outputFilename}.fasta`);
                    // console.log('output file removed');
                }
            );

        });*/

    }
}