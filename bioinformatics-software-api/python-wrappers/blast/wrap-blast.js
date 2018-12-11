/**
 * Created by anuradhawick on 9/3/17.
 */
import PythonShell from 'python-shell';

export default class WrapBlast {
    options;

    constructor(program = 'blastn', database = 'nr', sequence = null, outputFormat = 'Text') {
        this.options = {
            mode: 'json',
            pythonPath: process.env.PYTHON_PATH || 'python',
            pythonOptions: ['-u'],
            scriptPath: `${__dirname}`,
            args: [program, database, sequence, outputFormat]
        };
    }

    exec() {
        const shell = new PythonShell(`blast.py`, this.options);

        return new Promise((resolve, reject) => {
            shell.end((err) => {
                if (err) reject("BLAST error\n" + err);
            });

            shell.on('message', (message) => {
                // const lines = message.split('\n');
                // const count = lines.length;
                //
                // for (let i = 0; i < count; i++) {
                //     let line = lines.shift();
                //
                //     if(line.match(/<PRE>/)) {
                //         break;
                //     }
                // }

                // resolve({step: 'blast', output: lines.join('\n')});
                resolve({step: 'blast', output: message.output, text: message.text});
            });
        });

    }
}