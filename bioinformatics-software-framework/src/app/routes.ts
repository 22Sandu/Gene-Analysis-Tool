/**
 * Created by anuradhawick on 1/23/18.
 */
import { TutorialPageComponent } from './pages/tutorial-page/tutorial-page.component';
import { EditorViewComponent } from './pages/editor-view/editor-view.component';
import { AdminPageComponent } from './pages/admin-page/admin-page.component';
import { HomePageComponent } from './pages/home-page/home-page.component';
import { LoginPageComponent } from './pages/login-page/login-page.component';
import { CreateComponentComponent } from './pages/admin-page/create-component/create-component.component';
import { ManageComponentsComponent } from './pages/admin-page/manage-components/manage-components.component';
import { AboutUsComponent } from './pages/about-us/about-us.component';
import { FormsPageComponent } from './pages/forms-page/forms-page.component';

export const routes = [
    {
      path: '',
      component: EditorViewComponent,
      children: [
        {
          path: '',
          redirectTo: 'editor',
          pathMatch: 'full'
        },
        {
          path: 'editor',
          component: HomePageComponent,
        },
        {
          path: 'tutorial',
          component: TutorialPageComponent
        },
        {
          path: 'form',
          component: FormsPageComponent
        },
        {path: "about-us", component: AboutUsComponent},
        {
          path: 'admin',
          component: AdminPageComponent,
          children: [
            {
              path: '',
              redirectTo: 'create',
              pathMatch: 'full'
            },
            {
              path: 'create',
              component: CreateComponentComponent
            },
            {
              path: 'manage',
              component: ManageComponentsComponent
            }
          ]
        }
      ]
    },
    {
      path: 'login',
      component: LoginPageComponent
    }
  ]
;
