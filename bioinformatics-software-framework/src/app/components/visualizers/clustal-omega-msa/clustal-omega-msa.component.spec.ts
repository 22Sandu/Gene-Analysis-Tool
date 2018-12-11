import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { ClustalOmegaMsaComponent } from './clustal-omega-msa.component';

describe('ClustalOmegaMsaComponent', () => {
  let component: ClustalOmegaMsaComponent;
  let fixture: ComponentFixture<ClustalOmegaMsaComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ ClustalOmegaMsaComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ClustalOmegaMsaComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
