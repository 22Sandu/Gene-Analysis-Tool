import { Component, Input, OnInit } from '@angular/core';
import * as _ from 'lodash';

declare const $;

@Component({
  selector: 'app-home-page-help-modal',
  templateUrl: './home-page-help-modal.component.html',
  styleUrls: ['./home-page-help-modal.component.css']
})
export class HomePageHelpModalComponent implements OnInit {

  private helpStep: number = 1;

  constructor() {
  }

  ngOnInit() {
  }

  openHelpModal() {
    this.helpStep = 1;
    $("[modal='home-page-help']").modal('show');
  }

  prevHelpModalStep() {
    const stepVal = this.helpStep - 1;

    this.helpStep = _.max([stepVal, 1]);
  }

  nextHelpModalStep() {
    const stepVal = this.helpStep + 1;

    this.helpStep = _.min([stepVal, 3]);
  }

}
