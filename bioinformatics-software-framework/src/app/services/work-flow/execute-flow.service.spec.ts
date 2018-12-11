import { TestBed, inject } from '@angular/core/testing';

import { ExecuteFlowService } from './execute-flow.service';

describe('ExecuteFlowService', () => {
  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [ExecuteFlowService]
    });
  });

  it('should be created', inject([ExecuteFlowService], (service: ExecuteFlowService) => {
    expect(service).toBeTruthy();
  }));
});
