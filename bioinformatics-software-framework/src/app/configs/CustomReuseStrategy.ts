/**
 * Created by anuradhawick on 1/23/18.
 */
import { ActivatedRouteSnapshot, RouteReuseStrategy, DetachedRouteHandle } from '@angular/router';

interface RouteStorageObject {
  snapshot: ActivatedRouteSnapshot;
  handle: DetachedRouteHandle;
}

export class CustomReuseStrategy implements RouteReuseStrategy {

  storedRoutes: { [key: string]: RouteStorageObject } = {};

  shouldDetach(route: ActivatedRouteSnapshot): boolean {
    let detach: boolean = true;

    return detach;
  }

  store(route: ActivatedRouteSnapshot, handle: DetachedRouteHandle): void {
    let storedRoute: RouteStorageObject = {
      snapshot: route,
      handle: handle
    };

    this.storedRoutes[route.routeConfig.path] = storedRoute;
  }

  shouldAttach(route: ActivatedRouteSnapshot): boolean {

    let canAttach: boolean = !!route.routeConfig && !!this.storedRoutes[route.routeConfig.path];

    if (canAttach) {
      let willAttach: boolean = true;

      let paramsMatch: boolean = this.compareObjects(route.params, this.storedRoutes[route.routeConfig.path].snapshot.params);
      let queryParamsMatch: boolean = this.compareObjects(route.queryParams, this.storedRoutes[route.routeConfig.path].snapshot.queryParams);

      return paramsMatch && queryParamsMatch;
    } else {
      return false;
    }
  }

  retrieve(route: ActivatedRouteSnapshot): DetachedRouteHandle {
    if (!route.routeConfig || !this.storedRoutes[route.routeConfig.path]) return null;

    return this.storedRoutes[route.routeConfig.path].handle;
  }

  shouldReuseRoute(future: ActivatedRouteSnapshot, curr: ActivatedRouteSnapshot): boolean {
    return future.routeConfig === curr.routeConfig;
  }

  private compareObjects(base: any, compare: any): boolean {
    for (let baseProperty in base) {
      if (compare.hasOwnProperty(baseProperty)) {
        switch (typeof base[baseProperty]) {
          case 'object':
            if (typeof compare[baseProperty] !== 'object' || !this.compareObjects(base[baseProperty], compare[baseProperty])) {
              return false;
            }
            break;
          case 'function':
            if (typeof compare[baseProperty] !== 'function' || base[baseProperty].toString() !== compare[baseProperty].toString()) {
              return false;
            }
            break;
          default:
            if (base[baseProperty] != compare[baseProperty]) {
              return false;
            }
        }
      } else {
        return false;
      }
    }

    return true;
  }
}
