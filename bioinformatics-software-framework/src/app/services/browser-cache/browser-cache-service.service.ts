import { Injectable } from '@angular/core';
import * as _ from 'lodash';

@Injectable()
export class BrowserCacheService {

  constructor() {
  }

  saveComponentTree(tree: any) {
    localStorage.setItem('tree', JSON.stringify(tree));
  }

  getComponentTree() {
    let tree = localStorage.getItem('tree');

    return tree === null ? [] : JSON.parse(tree);
  }

}
