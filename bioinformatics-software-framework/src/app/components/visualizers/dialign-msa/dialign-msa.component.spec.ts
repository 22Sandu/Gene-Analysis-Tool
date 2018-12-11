import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { DialignMsaComponent } from './dialign-msa.component';

describe('DialignMsaComponent', () => {
  let component: DialignMsaComponent;
  let fixture: ComponentFixture<DialignMsaComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ DialignMsaComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(DialignMsaComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
