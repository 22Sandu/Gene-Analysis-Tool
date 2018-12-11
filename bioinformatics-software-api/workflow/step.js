/**
 * Created by anuradhawick on 9/11/17.
 */
export default class Step {
    nextStep = null;
    isInitial;
    task;
    data;

    constructor(task, data = null, isInitial = false) {
        this.isInitial = isInitial;
        this.task = task;
        this.data = data;
    }

    setNext(step) {
        this.nextStep = step;
    }

    async start(data = null) {
        const result = await this.task.execute(data || this.data);
        if (this.nextStep) {
            return await this.nextStep.start(result);
        } else {
            return result;
        }
    }
}