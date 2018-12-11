/**
 * Created by Vijini on 10/21/17.
 */

/**
 * This is the command line wrapper for the multiple alignment program DIALIGN2-2.
 * Can be downloaded from http://bibiserv.techfak.uni-bielefeld.de/dialign/welcome.html by going to Download tab
 */

// TODO Should setup DIALIGN2_DIR environment variable as given
// TODO Windows : set DIALIGN2_DIR=<parent_path>/dialign2_dir
// TODO OSX : export DIALIGN2_DIR=<parent_path>/dialign2_dir

import cmd from 'node-cmd';
import fs from 'fs';
import uuidv4 from 'uuid/v4';

export default class WrapDialign {
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

            // Format DIALIGN command to execute
            const command = `${__dirname}/sources/dialign2-2 -fa -fn ${__dirname}/${outputFilename} ${__dirname}/${inputFilename}.fasta`;

            // Execute DIALIGN command
            cmd.get(command, (err, data, stderr) => {
                if (err || stderr) {
                    reject("DIALIGN Error\n" + err);
                } else {
                    // Read .fasta result
                    const content = fs.readFileSync(`${__dirname}/${outputFilename}`, 'utf8');

                    resolve({step: 'dialign', output: content});
                }

                // Remove input file
                fs.existsSync(`${__dirname}/${inputFilename}.fasta`) && fs.unlinkSync(`${__dirname}/${inputFilename}.fasta`);

                // Remove output file
                fs.existsSync(`${__dirname}/${outputFilename}`) && fs.unlinkSync(`${__dirname}/${outputFilename}`);

                // Remove .fa file formed to format input file
                fs.existsSync(`${__dirname}/${outputFilename}.fa`) && fs.unlinkSync(`${__dirname}/${outputFilename}.fa`);
            });

        });

    }
}