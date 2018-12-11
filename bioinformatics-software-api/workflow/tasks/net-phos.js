/**
 * Created by anuradhawick on 9/11/17.
 */
import WrapNetPhos from '../../python-wrappers/net-phos/net-phos'

export default class NetPhos {

    constructor() {
    }

    execute(...taskParams) {
        const wrap = new WrapNetPhos(taskParams[0].output);

        return wrap.exec();
    }
}