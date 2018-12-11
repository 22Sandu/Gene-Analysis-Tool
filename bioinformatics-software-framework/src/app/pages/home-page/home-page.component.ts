import { Component, OnInit, ViewChild } from '@angular/core';

import { TreeViewComponent } from '../../components/tree-view/tree-view.component';
import { DrawingBoardComponent } from '../../components/drawing-board/drawing-board.component';
import { ExecuteFlowService } from '../../services/work-flow/execute-flow.service';
import { PairwiseBlastComponent } from '../../components/visualizers/pairwise-blast/pairwise-blast.component';
import { ClustalOmegaMsaComponent } from '../../components/visualizers/clustal-omega-msa/clustal-omega-msa.component';
import { DialignMsaComponent } from '../../components/visualizers/dialign-msa/dialign-msa.component';
import { TCoffeeMsaComponent } from '../../components/visualizers/t-coffee-msa/t-coffee-msa.component';
import { BarGraphComponent } from '../../components/visualizers/bar-graph/bar-graph.component';

import * as _ from 'lodash';

@Component({
  selector: 'app-home-page',
  templateUrl: './home-page.component.html',
  styleUrls: ['./home-page.component.css'],
  providers: [
    TreeViewComponent,
    DrawingBoardComponent,
    ExecuteFlowService,
    PairwiseBlastComponent
  ]
})
export class HomePageComponent implements OnInit {
  @ViewChild(DrawingBoardComponent) drawingBoard: DrawingBoardComponent;
  @ViewChild(PairwiseBlastComponent) pairwiseBlasterView: PairwiseBlastComponent;
  @ViewChild(ClustalOmegaMsaComponent) clustalOmegaView: ClustalOmegaMsaComponent;
  @ViewChild(DialignMsaComponent) dialignView: DialignMsaComponent;
  @ViewChild(TCoffeeMsaComponent) tCoffeeView: TCoffeeMsaComponent;
  @ViewChild(BarGraphComponent) barGraphView: BarGraphComponent;


  private activeTreeItem: any = null;
  private resultData: any = null;

  constructor(private executionService: ExecuteFlowService) {
  }

  ngOnInit() {

  }

  selectTree(item) {
    this.activeTreeItem = item;
  }

  addStep() {
    this.activeTreeItem.isStep && this.drawingBoard.addStep(this.activeTreeItem);
  }

  execute() {
    this.clearUI();

    const steps = this.drawingBoard.getStepSequence();

    this.executionService.executeFlow(steps).then((results) => {
      this.resultData = results;
      console.log('results came ', results);
      _.each(results, (result: any) => {
        console.log('name', result.name);
        if (result.name === 'Visualize Output') {
          console.log(result.computedValue);
          _.each(result.computedValue.data, (visualizeData:any) => {
            switch (visualizeData.step) {
              case 'blast':
                this.pairwiseBlasterView.render(visualizeData.text);
                break;
              case 'clustal-omega':
                console.log('data1',visualizeData.output);
                //this.barGraphView.render(visualizeData.output);
                this.clustalOmegaView.render(visualizeData.output);
                break;
              case 'clustal-omega-max-align':
                this.clustalOmegaView.render(visualizeData.output);
                break;
              case 'dialign':
                this.dialignView.render(visualizeData.output);
                break;
              case 't-coffee':
                this.tCoffeeView.render(visualizeData.output);
                break;
            }
          });
        }
      });
    });
  }

  clearUI() {
    this.pairwiseBlasterView.clear();
    this.clustalOmegaView.clear();
    this.dialignView.clear();
    this.tCoffeeView.clear();
    //this.barGraphView.clear();
  }

}
