import { Component, OnInit, ViewChild } from '@angular/core';
import { ModalComponent } from '../../components/modals/modal/modal.component';

@Component({
  selector: 'app-admin-page',
  templateUrl: './admin-page.component.html',
  styleUrls: ['./admin-page.component.css']
})
export class AdminPageComponent implements OnInit {
  @ViewChild('helloModal') modal: ModalComponent;

  constructor() {
  }

  ngOnInit() {
    setTimeout(() => {
      this.modal.openModal();
    }, 2000);
  }

}
