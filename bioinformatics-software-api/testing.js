/**
 * Created by anuradhawick on 1/13/18.
 */

import fs from 'fs';
import WrapTCoffee from './python-wrappers/t-coffee/wrap-t-coffee';
import WrapClustalOmega from './python-wrappers/clustal-omega/wrap-clustal-omega';
import WrapDialign from './python-wrappers/dialign/wrap-dialign';
import WrapBlast from './python-wrappers/blast/wrap-blast';

let startTime = Date.now();

let strMSA = [];
let testData = [
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000100.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000200.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000300.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000400.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000500.txt',
];

let testDataMSA = [
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000100.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000200.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000300.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000400.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000500.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000600.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000700.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000800.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00000900.txt',
    '/Users/anuradhawick/Desktop/blastp-queries/SequenceLength_00001000.txt',
];


// Generate MSA strings
for (let i = 0; i < 6; i++) {
    let tempStr = ``;
    for (let j = i; j < i + 5; j++) {
        tempStr += fs.readFileSync(testDataMSA[j]);
    }
    strMSA.push(tempStr);
}

let f = fs.readFileSync(testData[4]);

async function doBlast(seq) {
    const bls = new WrapBlast('blastp','env_nr',seq);
    return (await bls.exec()).output;
}

async function doDialign(seq) {
    const dal = new WrapDialign(seq);
    return await dal.exec();
}

async function doClustalOmega(seq) {
    const comega = new WrapClustalOmega(seq);
    return (await comega.exec()).output;
}

async function doTCoffee(seq) {
    const tc = new WrapTCoffee(seq);
    return await tc.exec();
}


// console.log(f.toString())
// console.log(strMSA[4]);

async function mainParallel(i) {
    startTime = Date.now();
    let x = doClustalOmega(strMSA[i]);
    let y = doTCoffee(strMSA[i]);
    let z = doDialign(strMSA[i]);
    await x;
    await y;
    await z;
    console.log(`parallel time is ${Date.now() - startTime}`);
}

async function mainSerial(i) {
    startTime = Date.now();
    let x = await doClustalOmega(strMSA[i]);
    console.log(`clustal serial time is ${Date.now() - startTime}`);
    startTime = Date.now();
    let y = await doTCoffee(strMSA[i]);
    console.log(`T-COFFEE serial time is ${Date.now() - startTime}`);
    startTime = Date.now();
    let z = await doDialign(strMSA[i]);
    console.log(`DIALIGN serial time is ${Date.now() - startTime}`);
}

async function msaBlsSerial(i) {
    startTime = Date.now();
    let x = await doClustalOmega(strMSA[i]);
    let y = await doBlast(fs.readFileSync(testData[i]));

    console.log(`serial time is ${Date.now() - startTime}`);
}

async function msaBlsParallel(i) {
    startTime = Date.now();
    let x = doClustalOmega(strMSA[i]);
    let y = doBlast(fs.readFileSync(testData[i]));

    await x;
    await y;
    console.log(`parallel time is ${Date.now() - startTime}`);
}
import os from 'os';

let initMemUsage = os.totalmem() - os.freemem();
console.log(`Initial memory - ${initMemUsage/1024/1024}`)

async function clsService(concurrency) {
    let startTime = Date.now();
    let ps = [];

    for (let i = 0; i < concurrency; i++) {
        ps.push(doDialign(strMSA[5]));
    }

    let additionaUsage = os.totalmem() - os.freemem();
    let doneDone = false;

    let int = setInterval(() => {
        if (doneDone) {
            clearInterval(int)
        }
        if (os.totalmem() - os.freemem() > additionaUsage) {
            additionaUsage = os.totalmem() - os.freemem();
        }
    }, 100);

    for (let i = 0; i < concurrency; i++) {
        await ps[i];
    }
    doneDone = true;

    console.log(`time is ${Date.now() - startTime} Memory used ${additionaUsage/1024/1024 - initMemUsage/1024/1024}`);
}


async function main() {
    // for (let i = 0; i < 6; i++) {
    //     console.log(`running tests for seq ${(i + 1) * 100}-${(i + 5) * 100}`)
    //     await mainParallel(i);
    //     await mainSerial(i);
    //     console.log('\n');
    // }

    for (let i = 1; i <= 10; i++) {
            // console.log(`running tests for seq ${(i + 1) * 100}-${(i + 5) * 100}`);
            console.log(`running concurrency for ${i}`);
            // await mainSerial(i);
            await clsService(i);
            // await msaBlsParallel(i);
            console.log('\n');
        }
    // let x  = await doBlast(fs.readFileSync(testData[0]));
    // console.log(x);
}

main();

// console.log(strMSA[5])
// require('dotenv').config();
// import TreeViewDBHandler from './db/treeview-db-handler';
//
// const obj = {
//     text: "User Input",
//     id: "I1",
//     precedence: "1",
//     isStep: true,
//     description: "Enter the sequence as a manual user input",
//     inputs: [
//         {
//             name: "Sequence",
//             type: "text",
//             value: ""
//         }
//     ]
// };
//
//
// TreeViewDBHandler.getAllViews((res) => {
//     console.log(res);
// });