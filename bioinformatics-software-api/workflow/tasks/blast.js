/**
 * Created by anuradhawick on 9/11/17.
 */
import WrapBlast from '../../python-wrappers/blast/wrap-blast';

export default class Blast {
    program = 'blastn';
    database = 'nr';

    constructor(program) {
        // this.program = program;
    }

    execute(...taskParams) {
        const wrap = new WrapBlast(this.program, this.database, taskParams[0].output);

        return wrap.exec();
    }
}