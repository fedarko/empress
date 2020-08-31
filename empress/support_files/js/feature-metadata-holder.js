define(["underscore", "util"], function (_, util) {
    /**
     * @class FeatureMetadataHolder
     *
     * @param {Array} featureMetadataColumns Columns of the feature metadata.
     *                Note: The order of this array should match the order of
     *                      the arrays which are the values of tipMetadata and
     *                      intMetadata. If no feature metadata was provided
     *                      when generating an Empress visualization, this
     *                      parameter should be [] (and recurringValues should
     *                      be [], and tipMetadata and intMetadata should be
     *                      {}s).
     * @param {Array} recurringValues Array of values that occur more than once
     *                                in the feature metadata (considering both
     *                                the tip and internal node metadata).
     * @param {Object} tipMetadata Feature metadata for tips in the tree.
     *                             The keys of this Object should be tips'
     *                             postorder positions in the tree, and the
     *                             values should be an array of feature
     *                             metadata values (of the same length as
     *                             featureMetadataColumns).
     * @param {Object} intMetadata Feature metadata for internal nodes in tree,
     *                             formatted analogously to tipMetadata. Since
     *                             this structure (as with tipMetadata) is
     *                             indexed by nodes' postorder positions in the
     *                             tree, metadata duplicated across internal
     *                             nodes with the same name should literally be
     *                             present multiple times in this structure.
     */
    function FeatureMetadataHolder(
        featureMetadataColumns,
        recurringValues,
        tipMetadata,
        intMetadata
    ) {
        /**
         * @type{Array}
         * Feature metadata column names.
         * @private
         */
        this._featureMetadataColumns = featureMetadataColumns;

        /**
         * @type{Array}
         * Values that occur more than once in the tip / internal node feature
         * metadata.
         * @private
         */
        this._recurringValues = recurringValues;

        /**
         * @type{Number}
         * Number of unique values in this._recurringValues.
         * @private
         */
        this._numRecurringValues = this._recurringValues.length;

        /**
         * @type{Object}
         * Feature metadata: keys are tree node postorder positions, and
         * values are arrays of length equal to featureMetadataColumns.length.
         * To make coloring the tree (in this.getUniqueInfo()) easier, we split
         * this up into tip and internal node feature metadata objects.
         * @private
         */
        this._tipMetadata = tipMetadata;
        this._intMetadata = intMetadata;
    }

    /**
     * Returns a copy of an array of the feature metadata columns.
     *
     * @return {Array}
     */
    FeatureMetadataHolder.prototype.getCols = function () {
        return _.clone(this._featureMetadataColumns);
    };

    /**
     * Returns the 0-based index of a column name in the f.m. columns.
     *
     * @param {String} cat Feature metadata column name.
     *
     * @return {Number}
     *
     * @throw {Error} If cat is not present in this._featureMetadataColumns.
     */
    FeatureMetadataHolder.prototype.getColIdx = function (cat) {
        var idx = _.indexOf(this._featureMetadataColumns, cat);
        if (idx < 0) {
            throw new Error(
                "Feature metadata column " + cat + " not present in data."
            );
        }
        return idx;
    };

    /**
     * Returns true if a node has tip feature metadata, false otherwise.
     *
     * (There isn't really anything stopping us from making another one of
     * these functions that checks the internal node metadata for a node.
     * None of the code needs to do that right now, though.)
     *
     * @param {Number} node Postorder position of a node in the tree
     *
     * @return {Boolean}
     */
    FeatureMetadataHolder.prototype.hasTipMetadata = function (node) {
        return _.has(this._tipMetadata, node);
    };

    /**
     * Returns the metadata value of a given node at a given column index.
     *
     * @param {Number} node Postorder position of a node in the tree
     * @param {Number} colIdx 0-indexed position of a feature metadata column
     *                        in this._featureMetadataColumns; this should have
     *                        been computed by this.getColIdx()
     *                        TODO refactor to use col name
     * @param {String} nodeType Should be one of "tip" or "int": "tip"
     *                          indicates that this is a tip node, and "int"
     *                          indicates that this is an internal node
     *
     * @return {String} Feature metadata value
     *
     * @throw {Error} If nodeType is not "tip" or "int", or if the node and/or
     *                column index are not present in the tip or internal node
     *                feature metadata
     */
    FeatureMetadataHolder.prototype.getValue = function (
        node,
        colIdx,
        nodeType
    ) {
        var initialVal;
        if (nodeType === "tip") {
            initialVal = this._tipMetadata[node][colIdx];
        } else if (nodeType === "int") {
            initialVal = this._intMetadata[node][colIdx];
        } else {
            throw new Error(
                "Unrecognized nodeType " + nodeType + " specified."
            );
        }
        if (_.isUndefined(initialVal)) {
            throw new Error(
                "Node " +
                    node +
                    " and/or col index " +
                    colIdx +
                    " not in " +
                    nodeType +
                    " feature metadata."
            );
        }
        return this._uncompressValue(initialVal);
    };

    /**
     * Returns all feature metadata values for a given node, in the same order
     * as this._featureMetadataColumns.
     *
     * Unlike getValue(), this doesn't raise an error if the node isn't present
     * in the feature metadata -- instead, it'll just return null in that case.
     *
     * @param {Number} node Postorder position of a node in the tree
     * @param {String} nodeType Should be one of "tip" or "int": "tip"
     *                          indicates that this is a tip node, and "int"
     *                          indicates that this is an internal node
     *
     * @return {String or null} Row of feature metadata values (uncompressed,
     *                          so no need to worry about recurring values
     *                          being replaced with numbers) if the node is
     *                          present in the feature metadata based on the
     *                          specified nodeType; null if the node is not
     *                          present in the feature metadata based on
     *                          nodeType
     *
     * @throw {Error} If nodeType is not "tip" or "int"
     */
    FeatureMetadataHolder.prototype.getRow = function (node, nodeType) {
        var row;
        if (nodeType === "tip") {
            row = this._tipMetadata[node];
        } else if (nodeType === "int") {
            row = this._intMetadata[node];
        } else {
            throw new Error(
                "Unrecognized nodeType " + nodeType + " specified."
            );
        }
        // If the node isn't in the specified f.m. type, return null
        if (_.isUndefined(row)) {
            return null;
        }
        return _.map(row, this._uncompressValue);
    };

    /**
     * Attempts to un-compress a value in the feature metadata, if needed.
     *
     * If a value occurs more than once in the metadata, it will have
     * been replaced in the metadata with a Number pointing to an index in
     * this._recurringValues. Since all other values should be stored as
     * Strings in the metadata, this uses the type of the value to determine
     * what to do.
     *
     * @param {String or Number} val A value in the feature metadata.
     *
     * @return {String} The input value, if it was of type String; or the
     *                  val-th value in this._recurringValues, if it was of
     *                  type Number.
     *
     * @throw {Error} If val is a Number but it isn't in an integer in the
     *                range [0, this._numRecurringValues).
     */
    FeatureMetadataHolder.prototype._uncompressValue = function (val) {
        if (typeof val === "number") {
            if (val >= this._numRecurringValues || val < 0) {
                throw new Error(
                    "Invalid recurring value compression: numerical " +
                        "value is " +
                        val +
                        ", which is out of the range " +
                        "[0, " +
                        this._numRecurringValues +
                        ")."
                );
            } else if (!Number.isInteger(val)) {
                throw new Error(
                    "Invalid recurring value compression: numerical " +
                        "value " +
                        val +
                        " is not an integer."
                );
            }
            return this._recurringValues[val];
        }
        return val;
    };

    /**
     * Retrieve unique value information for a feature metadata field.
     *
     * @param {String} cat The feature metadata column to find information for.
     *                     Must be present in this._featureMetadataColumns or
     *                     an error will be thrown.
     * @param {String} method Defines what feature metadata to check.
     *                        If this is "tip", then only tip-level feature
     *                        metadata will be used. If this is "all", then
     *                        this will use both tip and internal node feature
     *                        metadata. If this is anything else, this will
     *                        throw an error.
     * @return {Object} An object with two keys:
     *                  -sortedUniqueValues: maps to an Array of the unique
     *                   values in this feature metadata field, sorted using
     *                   util.naturalSort().
     *                  -uniqueValueToFeatures: maps to an Object which maps
     *                   the unique values in this feature metadata column to
     *                   an array of the node name(s) with each value.
     */
    FeatureMetadataHolder.prototype.getUniqueInfo = function (cat, method) {
        var scope = this;
        // In order to access feature metadata for a given node, we need to
        // find the 0-based index in this._featureMetadataColumns that the
        // specified f.m. column corresponds to.
        var fmIdx = this.getColIdx(cat);

        // The coloring method influences how much of the feature metadata
        // we'll look at. (While we're at it, validate the coloring method.)
        var fmObjs;
        if (method === "tip") {
            fmObjs = [this._tipMetadata];
        } else if (method === "all") {
            fmObjs = [this._tipMetadata, this._intMetadata];
        } else {
            throw 'F. metadata coloring method "' + method + '" unrecognized.';
        }
        // Produce a mapping of unique values in this feature metadata
        // column to an array of the node name(s) with each value.
        var uniqueValueToFeatures = {};
        _.each(fmObjs, function (mObj) {
            _.mapObject(mObj, function (fmRow, node) {
                // need to convert to integer
                node = parseInt(node);
                // This is loosely based on how BIOMTable.getObsBy() works.
                var fmVal = scope._uncompressValue(fmRow[fmIdx]);
                if (_.has(uniqueValueToFeatures, fmVal)) {
                    uniqueValueToFeatures[fmVal].push(node);
                } else {
                    uniqueValueToFeatures[fmVal] = [node];
                }
            });
        });

        var sortedUniqueValues = util.naturalSort(
            Object.keys(uniqueValueToFeatures)
        );
        return {
            sortedUniqueValues: sortedUniqueValues,
            uniqueValueToFeatures: uniqueValueToFeatures,
        };
    };

    return FeatureMetadataHolder;
});
