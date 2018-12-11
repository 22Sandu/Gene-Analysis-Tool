import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { PairwiseBlastComponent } from './pairwise-blast.component';

describe('PairwiseBlastComponent', () => {
  let component: PairwiseBlastComponent;
  let fixture: ComponentFixture<PairwiseBlastComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ PairwiseBlastComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(PairwiseBlastComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
