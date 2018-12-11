/**
 * Created by Vijini on 10/20/17.
 */

import cmd from 'node-cmd';
import fs from 'fs';
import uuidv4 from 'uuid/v4';

export default class WrapMaxAlignPearl {
    input;

    // TODO these are the parameters
    constructor(inp = null) {
        this.input = inp;
    }

    exec() {

        // Get input file name
        const inputFilename = uuidv4();

        // Create input file with input data
        fs.writeFileSync(`${__dirname}/${inputFilename}.fasta`, this.input);

        // Get output file name
        const outputFilename = uuidv4();

        return new Promise((resolve, reject) => {
            const command = `perl ${__dirname}/maxalign.pl ${__dirname}/${inputFilename}.fasta`;

            // Execute clustal omega command
            cmd.get(command, (err, data, stderr) => {

                    if (err || stderr) {
                        reject("Max Align Error\n" + err);
                    } else {
                        // Read .fasta result
                        const content = fs.readFileSync(`./heuristic.fsa`, 'utf8');
                        // console.log('CONTENT: ', content);

                        resolve({step: 'clustal-omega-max-align', output: content});
                    }

                    // Remove input file
                    fs.existsSync(`${__dirname}/${inputFilename}.fasta`) && fs.unlinkSync(`${__dirname}/${inputFilename}.fasta`);

                    // Remove output file
                    fs.existsSync(`./heuristic.fsa`) && fs.unlinkSync(`./heuristic.fsa`);
                    fs.existsSync(`./heuristic_exclude_headers.txt`) && fs.unlinkSync(`./heuristic_exclude_headers.txt`);
                    fs.existsSync(`./heuristic_include_headers.txt`) && fs.unlinkSync(`./heuristic_include_headers.txt`);
                    // console.log('output file removed');
                }
            );

        });

    }
}