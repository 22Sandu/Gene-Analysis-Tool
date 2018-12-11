import { Component, Input, OnInit } from '@angular/core';
import blasterjs from 'biojs-vis-blasterjs';

@Component({
  selector: 'app-pairwise-blast',
  templateUrl: './pairwise-blast.component.html',
  styleUrls: ['./pairwise-blast.component.css']
})
export class PairwiseBlastComponent implements OnInit {
  // @Input() data: any;

  constructor() {
  }

  ngOnInit() {

  }

  render(data) {
    const alignments = data;

    new blasterjs({
      string: alignments,
      multipleAlignments: "blast-multiple-alignments",
      alignmentsTable: "blast-alignments-table",
      singleAlignment: "blast-single-alignment"
    });
  }

  clear() {
    document.getElementById('blast-multiple-alignments').innerHTML = '';
    document.getElementById('blast-alignments-table').innerHTML = '';
    document.getElementById('blast-single-alignment').innerHTML = '';
  }
}
