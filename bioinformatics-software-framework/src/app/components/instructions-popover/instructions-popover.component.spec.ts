import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { InstructionsPopoverComponent } from './instructions-popover.component';

describe('InstructionsPopoverComponent', () => {
  let component: InstructionsPopoverComponent;
  let fixture: ComponentFixture<InstructionsPopoverComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ InstructionsPopoverComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(InstructionsPopoverComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
