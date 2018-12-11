/**
 * Created by anuradhawick on 9/11/17.
 */
import WrapSignalP from '../../python-wrappers/signal-p/signalp'

export default class SignalP {

    constructor() {
    }

    execute(...taskParams) {
        const wrap = new WrapSignalP(taskParams[0].output);

        return wrap.exec();
    }
}