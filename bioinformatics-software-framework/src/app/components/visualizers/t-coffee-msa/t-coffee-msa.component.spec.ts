import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { TCoffeeMsaComponent } from './t-coffee-msa.component';

describe('TCoffeeMsaComponent', () => {
  let component: TCoffeeMsaComponent;
  let fixture: ComponentFixture<TCoffeeMsaComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ TCoffeeMsaComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(TCoffeeMsaComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
