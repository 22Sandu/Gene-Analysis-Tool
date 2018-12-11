import { Injectable } from '@angular/core';
import { Http } from '@angular/http';

import 'rxjs/add/operator/toPromise';

import * as config from '../../configs/index';

@Injectable()
export class TreeViewService {

  constructor(private http: Http) {
  }

  getTree() {
    const headers: Headers = new Headers();

    headers.append('content-type', 'application/json');

    return this.http.get(
      `${config.apiUrl}/api/flow-tree`)
      .toPromise().then(response => response.json())
  }
}
