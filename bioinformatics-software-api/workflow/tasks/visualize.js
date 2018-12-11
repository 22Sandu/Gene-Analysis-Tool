/**
 * Created by anuradhawick on 9/11/17.
 */

export default class Visualize {

    execute(...taskParams) {
        return new Promise((resolve) => {
            const resultObject = {
                type: taskParams[0].step,
                command: 'visualize',
                data: taskParams[0]
            };

            resolve(resultObject);
        });
    }
}