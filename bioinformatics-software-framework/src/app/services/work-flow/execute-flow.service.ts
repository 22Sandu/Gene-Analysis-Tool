import { Injectable } from '@angular/core';
import { Http, Headers } from '@angular/http';

import 'rxjs/add/operator/toPromise';

import * as config from '../../configs/index';

@Injectable()
export class ExecuteFlowService {
  constructor(private http: Http) {
  }

  executeFlow(steps) {
    const headers: Headers = new Headers();

    headers.append('content-type', 'application/json');
    console.log('steps',steps);
    return this.http.post(
      `${config.apiUrl}/api/flow`,
      {
        steps: steps
      }, {
        headers: headers
      }).toPromise().then(response => response.json())
  }

}
