/**
 * Created by anuradhawick on 1/23/18.
 */
import express from "express";

import TreeViewDBHandler from '../db/treeview-db-handler';

const router = express.Router();

router.use((req, res, next) => {
    next()
});

router.get('/flow-tree', (req, res) => {
    res.header('content-type', 'application/json');

    TreeViewDBHandler.getAllViews((views) => {
        res.send(views);
    });
});


module.exports = router;