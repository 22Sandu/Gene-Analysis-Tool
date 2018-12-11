import { Component, OnInit, ViewChild } from '@angular/core';
import { ModalComponent } from '../../modals/modal/modal.component';

@Component({
  selector: 'app-dialign-msa',
  templateUrl: './dialign-msa.component.html',
  styleUrls: ['./dialign-msa.component.css']
})
export class DialignMsaComponent implements OnInit {
  @ViewChild('dialignViewRaw') modal: ModalComponent;
  data: string = null;
  active: boolean = false;

  constructor() { }

  ngOnInit() {
  }

  render(data) {
    this.data = data;
    this.active = true;
  }

  clear() {
    this.data = null;
    this.active = false;
  }

  viewRaw() {
    this.modal.openModal();
  }

}
