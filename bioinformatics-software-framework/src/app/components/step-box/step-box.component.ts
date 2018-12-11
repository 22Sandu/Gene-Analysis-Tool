import { Component, EventEmitter, Input, OnInit, Output } from '@angular/core';

declare const jsPlumb: any;
declare const $: any;

@Component({
  selector: 'app-step-box',
  templateUrl: './step-box.component.html',
  styleUrls: ['./step-box.component.css']
})
export class StepBoxComponent implements OnInit {
  @Input() step: any;
  @Output() boxClick: EventEmitter<any> = new EventEmitter();
  @Output() remove: EventEmitter<any> = new EventEmitter();

  constructor() {
  }

  ngOnInit() {
    // $(".item").resizable({
    //   resize : function(event, ui) {
    //     jsPlumb.repaint(ui.helper);
    //   },
    //   handles: "all"
    // });
  }

  onBoxClick(event) {
    event.stopPropagation();
    this.boxClick.emit(this.step);
  }

  onInputPress(event, id) {
    event.stopPropagation();
    $(`#${id}`).dialog();
  }

  onDialogCancel(input, id) {
    input.tempValue = input.value;
    $(`#${id}`).dialog('close');
  }

  onDialogOk(input, id) {
    input.value = input.tempValue;
    $(`#${id}`).dialog('close');
  }

  onStepRemove(event){
    event.stopPropagation();
    this.remove.emit(this.step);
  }

}
