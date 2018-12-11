/**
 * Created by anuradhawick on 2/26/18.
 */
import * as _ from 'lodash';

import * as Task from '../workflow/tasks';

export default function run(step) {
    switch (step.name) {
        case 'User Input':
            return (new Task.Data(step.params[0].value)).execute();
            break;
        case 'BLAST':
            // TODO get blast params
            return (new Task.Blast(`blast`)).execute(step.parentResults[0]);
            break;
        case 'Clustal Omega':
            return (new Task.ClustalOmega()).execute(step.parentResults[0]);
            break;
        case 'Visualize Output':
            return (new Task.Visualize()).execute(step.parentResults);
            break;
        case 'Max Align':
            return (new Task.MaxAling()).execute(step.parentResults[0]);
            break;
        case 'DIALIGN':
            return (new Task.Dialign()).execute(step.parentResults[0]);
            break;
        case 'T-Coffee':
            return (new Task.TCoffee()).execute(step.parentResults[0]);
            break;
        case 'NetPhos':
            return (new Task.NetPhos()).execute(step.parentResults[0]);
            break;
        case 'ProP':
            return (new Task.ProP()).execute(step.parentResults[0]);
            break;
        case 'SignalP':
            return (new Task.SignalP()).execute(step.parentResults[0]);
            break;
    }
}