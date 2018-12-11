import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { HomePageHelpModalComponent } from './home-page-help-modal.component';

describe('HomePageHelpModalComponent', () => {
  let component: HomePageHelpModalComponent;
  let fixture: ComponentFixture<HomePageHelpModalComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ HomePageHelpModalComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(HomePageHelpModalComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
