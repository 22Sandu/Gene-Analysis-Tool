/**
 * Created by anuradhawick on 9/11/17.
 */

export default class Data {
    data;

    constructor(data) {
        this.data = data;
    }

    execute() {
        return new Promise(resolve => resolve({step: 'blast', output: this.data}));
    }
}