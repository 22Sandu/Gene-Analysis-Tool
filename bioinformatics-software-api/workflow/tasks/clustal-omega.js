/**
 * Created by Vijini on 10/21/17.
 */
import WrapClustalOmega from '../../python-wrappers/clustal-omega/wrap-clustal-omega';

export default class ClustelOmega {

    execute(...taskParams) {
        const wrap = new WrapClustalOmega(taskParams[0].output);

        return wrap.exec();
    }
}