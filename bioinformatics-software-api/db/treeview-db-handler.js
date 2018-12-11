/**
 * Created by anuradhawick on 9/2/17.
 */
import {MongoClient}  from 'mongodb';

export default class TreeViewDBHandler {

    static insertView(viewObject, callback) {
        MongoClient.connect(process.env.MONGO_URL, (err, client) => {
            const db = client.db(process.env.MONGO_DB_NAME);
            const collection = db.collection('treeViewData');

            collection.insertOne(viewObject, (err, result) => {
                callback(result);
            });
            client.close();
        });
    }

    static getAllViews(callback) {
        MongoClient.connect(process.env.MONGO_URL, (err, client) => {
            const db = client.db(process.env.MONGO_DB_NAME);
            const collection = db.collection('treeViewData');

            collection.find({}).toArray((err, result) => {
                callback(result);
            });
            client.close();
        });
    }
}