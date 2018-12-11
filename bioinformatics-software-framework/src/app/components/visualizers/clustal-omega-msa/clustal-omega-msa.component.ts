import { Component, OnInit, ViewChild } from '@angular/core';
import { ModalComponent } from '../../modals/modal/modal.component';
import { Chart } from 'chart.js';

declare const msa: any;

@Component({
  selector: 'app-clustal-omega-msa',
  templateUrl: './clustal-omega-msa.component.html',
  styleUrls: ['./clustal-omega-msa.component.css']
})
export class ClustalOmegaMsaComponent implements OnInit {
  @ViewChild('clustalOmegaViewRaw') modal: ModalComponent;
  data: string = null;
  active: boolean = false;
  chart = [];

  constructor() {

  }

  ngOnInit() {
  }

  render(data) {
    /*const alignments = data;
    const seqs = msa.io.fasta.parse(alignments);
    const m = msa({
      el: document.getElementById('clustal-omega-msa'),
      seqs: seqs
    });

    this.data = data;
    m.render();
    this.active = true;*/
    console.log(data);
    this.chart = new Chart('clustal-omega-msa', {
      type: 'bar',
      data: {
        datasets: [
          {
            data: [10, 20, 30]
          }
        ]
      }
    });
  }

  clear() {
    this.active = false;
    document.getElementById('clustal-omega-msa').innerHTML = '';
  }

  viewRaw() {
    this.modal.openModal();
  }
}
