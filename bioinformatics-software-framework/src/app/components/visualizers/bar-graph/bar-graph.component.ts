import { Component, OnInit, ViewChild } from '@angular/core';
import { ModalComponent } from '../../modals/modal/modal.component';

declare const msa: any;

@Component({
  selector: 'app-bar-graph',
  templateUrl: './bar-graph.component.html',
  styleUrls: ['./bar-graph.component.css']
})
export class BarGraphComponent implements OnInit {
  @ViewChild('barGraphViewRaw') modal: ModalComponent;
  data: string = null;
  active: boolean = false;

  constructor() {

  }

  ngOnInit() {
  }

  render(data) {
    /*const alignments = data;
    const seqs = msa.io.fasta.parse(alignments);
    const m = msa({
      el: document.getElementById('bar-graph'),
      seqs: seqs
    });*/

    this.data = data;
    //m.render();
    this.active = true;
    //console.log(this.data);
  }

  clear() {
    this.active = false;
    document.getElementById('bar-graph').innerHTML = '';
  }

  viewRaw() {
    this.modal.openModal();
  }
}
