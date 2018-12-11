import { Component, OnInit } from '@angular/core';
import { StepBoxComponent } from '../step-box/step-box.component';
import * as _ from 'lodash';
import uuidv4 from 'uuid/v4';

declare const $: any;
declare const jsPlumb: any;

@Component({
  selector: 'app-drawing-board',
  templateUrl: './drawing-board.component.html',
  styleUrls: ['./drawing-board.component.css'],
  providers: [StepBoxComponent]
})
export class DrawingBoardComponent implements OnInit {
  private steps: Array<any> = [];
  private links: Array<any> = [];
  private connectorProps: any = {
    connector: ["Straight"],
    anchor: ["Left", "Right"],
    endpoint: "Dot"
  };

  constructor() {
  }

  ngOnInit() {
    // $(".drawing-surface").panzoom();
    jsPlumb.ready(() => {
      jsPlumb.importDefaults({
        Connector: ["Flowchart"]
      });

      jsPlumb.bind("connectionDetached", (info) => {
        const conn = info.connection;
        _.remove(this.links, (link) => {
          return link['connection'] === conn;
        });
      });
    })
  }

  updateUI() {
    // Check for DOM every 10ms and run plumber
    const int = setInterval(() => {
      let rendered = true;

      _.each(this.steps, (step) => {
        rendered = rendered && (document.getElementById(step.stepId) !== null);
      });

      if (rendered) {
        jsPlumb.ready(() => {
          // Updating boxes
          _.each(this.steps, (step) => {
            jsPlumb.draggable(step.stepId);
          });

          // Updating connections
          _.each(this.links, (link) => {
            // search if already exist
            if (_.isEmpty(link.connection)) {
              const connection = jsPlumb.connect({
                source: link.from.stepId,
                target: link.to.stepId,
                paintStyle: {strokeStyle: "grey", lineWidth: 3},
                endpointStyle: {fillStyle: "red", outlineColor: "black", radius: 5},
                overlays: [
                  ["Arrow", {width: 12, length: 12, location: 0.67}]
                ]
              }, this.connectorProps);
              link['connection'] = connection;
            }
          });
        });
        clearInterval(int);
      }
    }, 10);
  }

  addStep(step) {
    const localStep = _.cloneDeep(step);

    localStep['stepId'] = uuidv4();
    this.steps.push(localStep);
    this.updateUI();

  }

  clearClicks() {
    _.each(this.steps, (step) => {
      step['clicked'] = false;
    });
  }

  onBoxClick(event) {
    // if the click count is 2, clear all
    let clickSet = _.filter(this.steps, (step) => {
      return step.clicked;
    });

    if (clickSet.length === 2) {
      this.clearClicks();
    }
    // find step from list and update
    let step = _.find(this.steps, event);

    step['clicked'] = !step['clicked'];
    clickSet.push(step);

    // check if a link could be made
    if (clickSet.length === 2) {
      this.addLink(clickSet);
    }
  }

  onBoardClick() {
    this.clearClicks();
  }

  onStepRemove(step) {
    jsPlumb.ready(() => {
      jsPlumb.remove(step.stepId);
    });
    _.remove(this.steps, (localStep) => {
      return step.stepId === localStep.stepId;
    });
  }

  addLink(clickSet) {
    // sort by precedence of steps
    const sortedClicks = _.sortBy(clickSet, (click: any) => {
      return click.precedence;
    });

    if (_.findIndex(this.links, (link) => {
        return link.to.stepId === sortedClicks[1].stepId && link.from.stepId === sortedClicks[0].stepId
      }) === -1) {
      this.links.push({
        from: sortedClicks[0],
        to: sortedClicks[1]
      });
      this.updateUI();
    }
  }

  getStepSequence() {
    const steps = _.map(this.steps, (stepObj: any) => {
        return {
          stepId: stepObj.stepId,
          name: stepObj.text,
          params: _.map(stepObj.inputs, (paramObj: any) => {
            return {
              name: paramObj.name,
              value: paramObj.value
            }
          }),
          nextStepIds: [],
          parentStepIds: [],
          computedValue: null,
          parentResults: []
        }
    });

    _.each(this.links, (link) => {
      // get the parent step and child step
      const parentStep =_.find(steps, (step) => {
        return step.stepId === link.from.stepId
      });
      const childStep =_.find(steps, (step) => {
        return step.stepId === link.to.stepId
      });

      parentStep.nextStepIds.push(childStep.stepId);
      childStep.parentStepIds.push(parentStep.stepId);
    });

    return steps;
  }

  updateStepSequence(steps) {
    _.each(steps, (step) => {

    });
  }
}
