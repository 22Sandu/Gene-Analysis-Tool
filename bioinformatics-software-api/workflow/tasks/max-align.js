/**
 * Created by anuradhawick on 9/11/17.
 */
import WrapMaxAlignPearl from '../../python-wrappers/max-align/wrap-max-align-pearl';

export default class MaxAlign {

    constructor() {
    }

    execute(...taskParams) {
        const wrap = new WrapMaxAlignPearl(taskParams[0].output);

        return wrap.exec();
    }
}