/**
 * Created by anuradhawick on 10/14/17.
 */
import express from "express";
import * as _ from 'lodash';

import exec from '../workflow/step-executor';

import Step from '../workflow/step';
import * as Task from '../workflow/tasks';

const router = express.Router();

router.use((req, res, next) => {
    next()
});

router.post('/flow', async (req, res) => {
    res.header('content-type', 'application/json');

    const incompleteList = req.body.steps;
    const completedList = [];

    function hasRunnable() {
        const result = _.find(incompleteList, (step) => {
            return step.parentStepIds.length === 0 || step.parentStepIds.length === step.parentResults.length;
        });
        return !_.isEmpty(result);
    }

    function getRunnables() {
        return _.remove(incompleteList, (step) => {
            return step.parentStepIds.length === 0 || step.parentStepIds.length === step.parentResults.length;
        });
    }

    async function executeStep(step) {
        // execute and obtain result
        const result = await exec(step);

        step.computedValue = result;

        // propagate results for children
        const nextStepIds = step.nextStepIds;
        _.each(nextStepIds, (id) => {
            const stepObj = _.find(incompleteList, (obj) => id === obj.stepId);
            stepObj.parentResults.push(step.computedValue);
        });
    }

    while (hasRunnable()) {
        const runnableList = getRunnables();
        let step = null;

        for (let i = 0; i < runnableList.length; i++) {
            step = runnableList[i];
            // console.log('running step', step.name);
            await executeStep(step);

            console.log(step.name);
            completedList.push(step);
        }
    }
    res.send(_.map(completedList, (step) => {
        return {
            name: step.name,
            stepId: step.stepId,
            computedValue: step.computedValue
        }
    }));
    // console.log('completed ', completedList)
});


module.exports = router;