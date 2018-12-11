/**
 * Created by anuradhawick on 9/11/17.
 */
import WrapProP from '../../python-wrappers/prop/prop'

export default class ProP {

    constructor() {
    }

    execute(...taskParams) {
        const wrap = new WrapProP(taskParams[0].output);

        return wrap.exec();
    }
}