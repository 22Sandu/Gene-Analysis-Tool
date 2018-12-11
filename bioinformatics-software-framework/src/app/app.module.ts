import { BrowserModule } from '@angular/platform-browser';
import { NgModule, NO_ERRORS_SCHEMA } from '@angular/core';
import { RouterModule, RouteReuseStrategy }   from '@angular/router';

import { AppComponent } from './app.component';
import { HomePageComponent } from './pages/home-page/home-page.component';
import { LoginPageComponent } from './pages/login-page/login-page.component';
import { TreeViewComponent } from './components/tree-view/tree-view.component';
import { DrawingBoardComponent } from './components/drawing-board/drawing-board.component';
import { StepBoxComponent } from './components/step-box/step-box.component';
import { OptionsComponent } from './components/options/options.component';
import { FormsModule } from '@angular/forms';
import { HttpModule } from '@angular/http';
import { PairwiseBlastComponent } from './components/visualizers/pairwise-blast/pairwise-blast.component';
import { ClustalOmegaMsaComponent } from './components/visualizers/clustal-omega-msa/clustal-omega-msa.component';
import { BarGraphComponent } from './components/visualizers/bar-graph/bar-graph.component';
import { InstructionsPopoverComponent } from './components/instructions-popover/instructions-popover.component';
import { HomePageHelpModalComponent } from './components/modals/home-page-help-modal/home-page-help-modal.component';
import { TestDataModalComponent } from './components/modals/test-data-modal/test-data-modal.component';
import { ModalComponent } from './components/modals/modal/modal.component';
import { TutorialPageComponent } from './pages/tutorial-page/tutorial-page.component';
import { EditorViewComponent } from './pages/editor-view/editor-view.component';
import { AdminPageComponent } from './pages/admin-page/admin-page.component';
import { CustomReuseStrategy } from './configs/CustomReuseStrategy';
import {routes} from './routes';
import { CreateComponentComponent } from './pages/admin-page/create-component/create-component.component';
import { ManageComponentsComponent } from './pages/admin-page/manage-components/manage-components.component';
import { BrowserCacheService } from './services/browser-cache/browser-cache-service.service';
import { AboutUsComponent } from './pages/about-us/about-us.component';
import { DialignMsaComponent } from './components/visualizers/dialign-msa/dialign-msa.component';
import { TCoffeeMsaComponent } from './components/visualizers/t-coffee-msa/t-coffee-msa.component';
import { EscapeHtmlPipe } from './configs/EscapeHtmlPipe';
import { FormsPageComponent } from './pages/forms-page/forms-page.component';

@NgModule({
  declarations: [
    AppComponent,
    HomePageComponent,
    LoginPageComponent,
    TreeViewComponent,
    DrawingBoardComponent,
    StepBoxComponent,
    OptionsComponent,
    PairwiseBlastComponent,
    ClustalOmegaMsaComponent,
    BarGraphComponent,
    InstructionsPopoverComponent,
    HomePageHelpModalComponent,
    TestDataModalComponent,
    ModalComponent,
    TutorialPageComponent,
    EditorViewComponent,
    AdminPageComponent,
    CreateComponentComponent,
    ManageComponentsComponent,
    AboutUsComponent,
    DialignMsaComponent,
    TCoffeeMsaComponent,
    EscapeHtmlPipe,
    FormsPageComponent
  ],
  imports: [
    BrowserModule,
    RouterModule.forRoot(routes),
    FormsModule,
    HttpModule
  ],
  providers: [
    {provide: RouteReuseStrategy, useClass: CustomReuseStrategy},
    BrowserCacheService
  ],
  bootstrap: [AppComponent],
  schemas: [NO_ERRORS_SCHEMA]
})
export class AppModule {
}


