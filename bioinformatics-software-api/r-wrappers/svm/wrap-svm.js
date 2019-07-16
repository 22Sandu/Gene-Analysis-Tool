var R = require("r-script");

export default class WrapSVM{

    constructor(){}

    exec(){
        var out = R("D:/Academics/Research/JBCB/Tool/Gene-Analysis-Tool/bioinformatics-software-api/r-wrappers/svm/svm.R")
        .data()
        .callSync();
        console.log(out.toString());
  
    }
}