import { Component, OnInit, ViewChild } from '@angular/core';
import { ModalComponent } from '../../modals/modal/modal.component';

@Component({
  selector: 'app-t-coffee-msa',
  templateUrl: './t-coffee-msa.component.html',
  styleUrls: ['./t-coffee-msa.component.css']
})
export class TCoffeeMsaComponent implements OnInit {
  @ViewChild('tCoffeeViewRaw') modal: ModalComponent;
  data: any = null;
  active: boolean = false;

  constructor() {
  }

  ngOnInit() {
  }

  render(data) {
    this.data = data;
    this.active = true;
  }

  clear() {
    this.data = '';
    this.active = false;
  }

  viewRaw() {
    this.modal.openModal();
  }

}
