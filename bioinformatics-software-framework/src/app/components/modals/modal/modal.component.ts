import { Component, Input, OnInit } from '@angular/core';
import uuidv4 from 'uuid/v4';

declare const $;

@Component({
  selector: 'app-modal',
  templateUrl: './modal.component.html',
  styleUrls: ['./modal.component.css']
})
export class ModalComponent implements OnInit {
  @Input() width: string = '50%';
  modalId: string;

  constructor() {
    this.modalId = uuidv4();
  }

  ngOnInit() {
  }

  openModal() {
    $(`[modal='modal-${this.modalId}']`).modal('show');
  }

}
