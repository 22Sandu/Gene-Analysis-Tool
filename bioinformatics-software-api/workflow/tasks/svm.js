import WrapSVM from '../../r-wrappers/svm/wrap-svm';

export default class Blast {
    // program = 'blastn';
    // database = 'nr';

    constructor(program) {
        // this.program = program;
    }

    execute(...taskParams) {
        const wrap = new WrapSVM();

        return wrap.exec();
    }
}