/**
 * Created by anuradhawick on 13/02/18.
 */
import WrapDialign from "../../python-wrappers/dialign/wrap-dialign";

export default class Dialign {

    execute(...taskParams) {
        const wrap = new WrapDialign(taskParams[0]);

        return wrap.exec();
    }
}