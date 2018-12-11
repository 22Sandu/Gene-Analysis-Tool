/**
 * Created by anuradhawick on 13/02/18.
 */

import WrapTCoffee from "../../python-wrappers/t-coffee/wrap-t-coffee";

export default class TCoffee {

    execute(...taskParams) {
        const wrap = new WrapTCoffee(taskParams[0].output);

        return wrap.exec();
    }
}