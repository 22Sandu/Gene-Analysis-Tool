/**
 * Created by anuradhawick on 2/24/18.
 */

export class Step {
  inputs: Array<any> = [];
  output: any = null;
  noOfParents: number = 0;
  type: string = null;

  constructor() {

  }

  setInput(input) {
    this.inputs.push(input);
    this.noOfParents++;
  }
}

export class Graph {

  constructor() {

  }


}
