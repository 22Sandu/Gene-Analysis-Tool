import { Component, OnInit } from '@angular/core';

declare const $;

@Component({
  selector: 'app-test-data-modal',
  templateUrl: './test-data-modal.component.html',
  styleUrls: ['./test-data-modal.component.css']
})
export class TestDataModalComponent implements OnInit {

  constructor() { }

  ngOnInit() {
  }

  openHelpModal() {
    $("[modal='home-page-test-data']").modal('show');
  }

}
