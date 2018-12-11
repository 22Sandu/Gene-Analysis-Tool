import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { PopoverTooltipComponent } from './popover-tooltip.component';

describe('PopoverTooltipComponent', () => {
  let component: PopoverTooltipComponent;
  let fixture: ComponentFixture<PopoverTooltipComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ PopoverTooltipComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(PopoverTooltipComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
