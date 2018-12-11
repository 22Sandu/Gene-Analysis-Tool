import { Component, OnInit, ViewChild } from '@angular/core';
import { HomePageHelpModalComponent } from '../../components/modals/home-page-help-modal/home-page-help-modal.component';
import { TestDataModalComponent } from '../../components/modals/test-data-modal/test-data-modal.component';
import * as _ from  'lodash';

@Component({
  selector: 'app-editor-view',
  templateUrl: './editor-view.component.html',
  styleUrls: ['./editor-view.component.css']
})
export class EditorViewComponent implements OnInit {
  @ViewChild(HomePageHelpModalComponent) homePageHelpModal: HomePageHelpModalComponent;
  @ViewChild(TestDataModalComponent) testDataModal: TestDataModalComponent;

  // localStorage = localStorage;

  constructor() { }

  ngOnInit() {
    // if (this.getParameterByName("type") === 'tutorial' && _.isEmpty(localStorage.getItem("type"))) {
    //   localStorage.setItem("type", "tutorial");
    // } else if (_.isEmpty(localStorage.getItem("type"))) {
    //   localStorage.setItem("type", "tool");
    // }
  }

  openHelpModal() {
    this.homePageHelpModal.openHelpModal();
  }

  openTestDataModal() {
    this.testDataModal.openHelpModal();
  }

  getParameterByName(name: string) {
    const searchParams = new URLSearchParams(window.location.search);

    return searchParams.get(name) || null;
  }
}
