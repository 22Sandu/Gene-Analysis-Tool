import { TestBed, inject } from '@angular/core/testing';

import { BrowserCacheService } from './browser-cache-service.service';

describe('BrowserCacheService', () => {
  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [BrowserCacheService]
    });
  });

  it('should be created', inject([BrowserCacheService], (service: BrowserCacheService) => {
    expect(service).toBeTruthy();
  }));
});
