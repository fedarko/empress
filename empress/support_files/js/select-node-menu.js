define(["underscore", "util"], function (_, util) {
    function SelectedNodeMenu(empress, drawer) {
        this.empress = empress;
        this.drawer = drawer;

        // General elements
        this.box = document.getElementById("menu-box");
        this.nodeNameLabel = document.getElementById("menu-box-node-id");
        this.nodeNameWarning = document.getElementById(
            "menu-box-node-name-warning"
        );
        this.nodeNotInTableWarning = document.getElementById(
            "menu-box-node-not-in-table-warning"
        );
        this.nodeLengthContainer = document.getElementById(
            "menu-box-node-length-container"
        );
        this.nodeLengthLabel = document.getElementById("menu-box-node-length");

        // Sample metadata elements
        this.smSection = document.getElementById("menu-sm-section");
        this.smTable = document.getElementById("menu-sm-table");
        this.smHeader = document.getElementById("menu-sm-header");
        this.sel = document.getElementById("menu-select");
        this.addBtn = document.getElementById("menu-add-btn");
        this.smAddSection = document.getElementById("menu-sm-add-section");
        this.smNotes = document.getElementById("menu-box-notes");

        // Feature metadata elements
        this.fmTable = document.getElementById("menu-fm-table");
        this.fmHeader = document.getElementById("menu-fm-header");
        this.fmSection = document.getElementById("menu-fm-section");

        this.nodeKeys = null;
        this.fields = [];
        this.hiddenCallback = null;
        this.visibleCallback = null;
        this._samplesInSelection = [];
        this.initialize();
    }

    /**
     * Un-hides a HTMLElement.
     *
     * @param {HTMLElement} ele
     */
    function show(ele) {
        ele.classList.remove("hidden");
    }

    /**
     * Sets the textContent of a HTMLElement and un-hides it.
     *
     * @param {HTMLElement} warningEle
     * @param {String} msg
     */
    function updateAndShow(ele, msg) {
        ele.textContent = msg;
        show(ele);
    }

    /**
     * Hides a HTMLElement.
     *
     * @param {HTMLElement} ele
     */
    function hide(ele) {
        ele.classList.add("hidden");
    }

    /**
     * Initializes the state machine. Adds metadata field options to drop down
     * menu, and creates the add button click event.
     */
    SelectedNodeMenu.prototype.initialize = function () {
        var scope = this;

        if (this.empress.isCommunityPlot) {
            // add items to select
            var selOpts = this.empress.getSampleCategories();
            for (var i = 0; i < selOpts.length; i++) {
                var opt = document.createElement("option");
                opt.value = selOpts[i];
                opt.innerHTML = selOpts[i];
                this.sel.appendChild(opt);
            }
            // add event to add button
            var click = function () {
                var val = scope.sel.value;
                scope.sel.options[scope.sel.selectedIndex].remove();
                // Hide the add button and related elements when all fields
                // are added: https://github.com/biocore/empress/issues/272
                if (scope.sel.options.length === 0) {
                    hide(scope.smAddSection);
                }
                scope.fields.push(val);
                show(scope.smHeader);
                scope.showNodeMenu();
            };
            this.addBtn.onclick = click;
        } else {
            hide(this.smSection);
        }
    };

    /*
     * Creates a HTML table describing sample presence info for a feature.
     *
     * This is set up as a static method
     * (https://stackoverflow.com/a/1635143/10730311) to make testing easier
     * (and also because it really doesn't need to depend on the state of this
     * object).
     *
     * @param{Object} ctData Two-dimensional mapping: The keys are the
     *                       sample metadata fields to include in the table,
     *                       and the values are Objects mapping unique values
     *                       in these fields to numbers describing this
     *                       feature's presence for these values.
     *                       e.g. {"body-site": {"gut": 5, "tongue": 2}}
     */
    SelectedNodeMenu.prototype.makeSampleMetadataTable = function (ctData) {
        // loop over all metadata fields the user has decided to show
        var sortedFields = util.naturalSort(_.keys(ctData));
        for (var i = 0; i < sortedFields.length; i++) {
            var field = sortedFields[i];

            // Create new rows in menu-table: the first row is for this
            // metadata field's "headers" (the unique values in the field,
            // e.g. "gut", "tongue", etc. for a field like body site), and
            // the second row is for the sample presence data for
            // the selected tree node within these unique values.
            //
            // Each group of two rows additionally has a header cell
            // on its leftmost side which spans both the header and data
            // rows; this header cell contains the name of the selected
            // metadata field, and has some fancy CSS that keeps it frozen
            // in place as the user scrolls the table horizontally.
            var fieldHeaderRow = this.smTable.insertRow(-1);
            var fieldHeaderCell = fieldHeaderRow.insertCell(-1);
            fieldHeaderCell.innerHTML = "<strong>" + field + "</strong>";
            fieldHeaderCell.rowSpan = 2;
            fieldHeaderCell.classList.add("menu-box-header-cell");
            fieldHeaderCell.classList.add("frozen-cell");

            var fieldDataRow = this.smTable.insertRow(-1);

            // add row values for this metadata field, one column at a time
            var categories = util.naturalSort(_.keys(ctData[field]));
            for (var j = 0; j < categories.length; j++) {
                var categoryHeaderCell = fieldHeaderRow.insertCell(-1);
                categoryHeaderCell.innerHTML =
                    "<strong>" + categories[j] + "</strong>";
                var categoryDataCell = fieldDataRow.insertCell(-1);
                categoryDataCell.innerHTML = ctData[field][categories[j]];
            }
        }
    };

    /*
     * Creates a HTML table (and a header) describing feature metadata.
     *
     * This checks to make sure that there actually is feature metadata (and
     * that the requested node has feature metadata) before creating things --
     * unlike makeSampleMetadataTable(), it's expected that some nodes may not
     * have any feature metadata information, or that feature metadata may not
     * have even been provided in the first place. (If this is the case, this
     * function will hide the fmHeader and fmTable elements.)
     *
     * @param{String} nodeName Name of the node to create this table for.
     *                         Duplicate names (for internal nodes) are ok.
     * @param{String} tipOrInt "tip" to query tip metadata, "int" to query
     *                         internal node metadata. Other values will cause
     *                         an error.
     */
    SelectedNodeMenu.prototype.makeFeatureMetadataTable = function (
        nodeName,
        tipOrInt
    ) {
        var fmShown = false;
        var fmCols = this.empress.getFeatureMetadataCategories();
        if (fmCols.length > 0) {
            var mdObj;
            // TODO: it'd be nice to add methods to Empress that return the
            // feature metadata, so we can avoid having to access "private"
            // attributes here. This is planned as part of the
            // https://github.com/biocore/empress/issues/337 refactoring.
            if (tipOrInt === "tip") {
                mdObj = this.empress._tipMetadata;
            } else if (tipOrInt === "int") {
                mdObj = this.empress._intMetadata;
            } else {
                throw new Error("Invalid tipOrInt value: " + tipOrInt);
            }
            if (_.has(mdObj, nodeName)) {
                var headerRow = this.fmTable.insertRow(-1);
                var featureRow = this.fmTable.insertRow(-1);
                for (var x = 0; x < fmCols.length; x++) {
                    var colName = fmCols[x];
                    var colCell = headerRow.insertCell(-1);
                    colCell.innerHTML = "<strong>" + colName + "</strong>";
                    var dataCell = featureRow.insertCell(-1);
                    dataCell.innerHTML = mdObj[nodeName][x];
                }
                show(this.fmSection);
                fmShown = true;
            }
        }
        if (!fmShown) {
            hide(this.fmSection);
        }
    };

    /**
     * Displays the node selection menu. nodeKeys must be set in order to use
     * this method.
     */
    SelectedNodeMenu.prototype.showNodeMenu = function () {
        // make sure the state machine is set
        if (this.nodeKeys === null) {
            throw "showNodeMenu(): Nodes have not been selected.";
        }

        // grab the name of the node
        var emp = this.empress;
        var nodeKeys = this.nodeKeys;
        var node = nodeKeys[0];
        var name = emp.getNodeInfo(node, "name");

        if (name === null) {
            this.nodeNameLabel.textContent = "Unnamed node";
        } else {
            this.nodeNameLabel.textContent = "Name: " + name;
        }

        this.smNotes.textContent = "";
        hide(this.nodeNameWarning);
        hide(this.nodeNotInTableWarning);

        // show either leaf or internal node
        var t = emp._tree;
        if (t.isleaf(t.postorderselect(this.nodeKeys[0]))) {
            this.showLeafNode();
        } else {
            this.showInternalNode();
        }

        // place menu-node menu next to node
        // otherwise place the (aggregated) node-menu over the root of the tree
        this.updateMenuPosition();

        // show table
        show(this.box);

        if (this.visibleCallback !== null) {
            this.visibleCallback(this._samplesInSelection);
        }
    };

    /**
     * Creates the node menu-table for a tip node. nodeKeys must be set in
     * before this function is called.
     */
    SelectedNodeMenu.prototype.showLeafNode = function () {
        // test to make sure nodeKeys is set
        if (this.nodeKeys === null) {
            throw "showLeafNode(): nodeKeys is not set!";
        }

        // test to make sure the leaf node is unique
        // (This should already be enforced in the Python side of things, but
        // we may as well be extra cautious.)
        if (this.nodeKeys.length > 1) {
            throw "showLeafNode(): Leaf nodes must be unique!";
        }

        // get the name of the tip
        var node = this.nodeKeys[0];

        // 1. Add feature metadata information (if present for this tip; if
        // there isn't feature metadata for this tip, the f.m. UI elements in
        // the selected node menu will be hidden)
        this.makeFeatureMetadataTable(node, "tip");

        this.setNodeLengthLabel(node);

        // 2. Add sample presence information for this tip (only if this data
        // is available in the first place, and if the user has selected at
        // least one field to show sample presence information for)
        if (this.empress.isCommunityPlot) {
            var ctData = this.empress.computeTipSamplePresence(
                node,
                this.fields
            );
            if (_.isNull(ctData)) {
                // This tip isn't present in the table, so we don't have sample
                // presence information for it.
                updateAndShow(
                    this.nodeNotInTableWarning,
                    "This is a tip in the tree. However, it is not " +
                        "present in the input feature table, so we cannot " +
                        "show sample presence information for it."
                );
                hide(this.smSection);
            } else {
                // 2.1 The samples represented by this tip are sent to Emperor.

                // Check if this tip is present in the BIOM table. The array
                // returned by BIOMTable.getObsIDsDifference() contains the
                // feature IDs present in the input array but not in the BIOM
                // table -- so if the length of this array is zero, this
                // feature is present in the table.
                var diff = this.empress._biom.getObsIDsDifference([node]);
                if (diff.length == 0) {
                    this._samplesInSelection = this.empress._biom.getSamplesByObservations(
                        [node]
                    );
                } else {
                    this._samplesInSelection = [];
                }
                this._checkTips(diff);

                this.makeSampleMetadataTable(ctData);
                if (this.fields.length > 0) {
                    this.smNotes.textContent =
                        "This is a tip in the tree. These values " +
                        "represent the number of unique samples that " +
                        "contain this node.";
                }
                show(this.smSection);
            }
        }
    };

    /**
     * Creates the node menu-table for internal nodes. nodeKeys must be set in
     * before this function is called. Furthermore, if there are more than key
     * in nodeKeys, then the keys must represent internal nodes with the same
     * name in the newick tree.
     */
    SelectedNodeMenu.prototype.showInternalNode = function () {
        // test to make sure nodeKeys is set
        if (this.nodeKeys === null) {
            throw "showInternalNode(): nodeKeys is not set!";
        }

        var name = this.empress.getNodeInfo(this.nodeKeys[0], "name");

        // Figure out whether or not we know the actual node in the tree (for
        // example, if the user searched for a node with a duplicate name, then
        // we don't know which node the user was referring to). This impacts
        // whether or not we show the sample presence info for this node.
        var isUnambiguous = this.nodeKeys.length === 1;

        // This is not necessarily equal to this.nodeKeys. If an internal node
        // with a duplicate name was clicked on then this.nodeKeys will only
        // have a single entry (the node that was clicked on): but
        // keysOfNodesWithThisName will accordingly have multiple entries.
        // The reason we try to figure this out here is so that we can
        // determine whether or not to show a warning about duplicate names
        // in the menu.
        if (name !== null) {
            var keysOfNodesWithThisName = this.empress._tree.getNodesWithName(
                name
            );
            if (keysOfNodesWithThisName.length > 1) {
                updateAndShow(
                    this.nodeNameWarning,
                    "Warning: " +
                        keysOfNodesWithThisName.length +
                        " nodes exist with the above name."
                );
            }
        } else {
            updateAndShow(
                this.nodeNameWarning,
                "No name was provided for this node in the input tree file."
            );
        }

        // 1. Add feature metadata information (if present) for this node
        // (Note that we allow duplicate-name internal nodes to have
        // feature metadata; this isn't a problem)
        this.makeFeatureMetadataTable(this.nodeKeys[0], "int");

        // 2. Compute sample presence information for this node.
        // (NOTE: this does not prevent "double-counting" samples, so the
        // aggregation for duplicate names should be fixed.)

        // force-reset the selection buffer
        this._samplesInSelection = [];

        if (isUnambiguous) {
            // this.nodeKeys has a length of 1
            var nodeKey = this.nodeKeys[0];
            this.setNodeLengthLabel(nodeKey);
            if (this.empress.isCommunityPlot) {
                var tips = this.empress._tree.findTips(nodeKey);

                var emp = this.empress;
                var samplePresence = emp.computeIntSamplePresence(
                    nodeKey,
                    this.fields
                );

                if (_.isNull(samplePresence.fieldsMap)) {
                    updateAndShow(
                        this.nodeNotInTableWarning,
                        "This is an internal node in the tree. None of " +
                            "this node's descendant tips are present in the " +
                            "input feature table, so we cannot show sample " +
                            "presence information for it."
                    );
                    hide(this.smSection);
                } else {
                    // used for the emperor callback
                    this._samplesInSelection = this._samplesInSelection.concat(
                        samplePresence.samples
                    );
                    if (this.fields.length > 0) {
                        this.makeSampleMetadataTable(samplePresence.fieldsMap);
                        this.smNotes.textContent =
                            "This is an internal node in the tree. These " +
                            "values represent the number of unique samples that " +
                            "contain any of this node's descendant tips.";
                        show(this.smSection);
                    }
                }
                this._checkTips(samplePresence.diff);
            }
        } else {
            // If isUnambiguous is false, no notes will be shown and the sample
            // presence info (including the table and notes) will be hidden
            hide(this.smSection);
            hide(this.nodeLengthContainer);
        }
    };

    /**
     * Given an array of tip names that are not present in the BIOM table,
     * warns the user about them using a toast message.
     */
    SelectedNodeMenu.prototype._checkTips = function (diff) {
        if (
            diff.length &&
            (this.visibleCallback !== null || this.hiddenCallback !== null)
        ) {
            util.toastMsg(
                "The following tips are not represented by your " +
                    "feature table and ordination: " +
                    diff.join(", ")
            );
        }
    };

    /**
     * Updates and shows the node length UI elements for a given node.
     *
     * (If the node is the root of the tree, this will actually hide the UI
     * elements. See Empress.getNodeLength() for details.)
     *
     * @param {Number} nodeKey Postorder position of a node in the tree.
     */
    SelectedNodeMenu.prototype.setNodeLengthLabel = function (nodeKey) {
        var nodeLength = this.empress.getNodeLength(nodeKey);
        if (nodeLength !== null) {
            this.nodeLengthLabel.textContent = nodeLength;
            show(this.nodeLengthContainer);
        } else {
            // Don't show the length for the root node
            hide(this.nodeLengthContainer);
        }
    };

    /**
     * Resets the state machine.
     */
    SelectedNodeMenu.prototype.clearSelectedNode = function () {
        this.smTable.innerHTML = "";
        this.nodeKeys = null;
        hide(this.box);
        hide(this.fmSection);
        this.fmTable.innerHTML = "";
        this.drawer.loadSelectedNodeBuff([]);
        this.empress.drawTree();

        if (this.hiddenCallback !== null) {
            this.hiddenCallback(this._samplesInSelection);
        }
        this._samplesInSelection = [];
    };

    /**
     * Sets the nodeKeys parameter of the state machine. This method will also
     * set the buffer to highlight the selected nodes.
     *
     * @param {Array} nodeKeys An array of node keys representing the
     *                         nodes to be selected. The keys should be the
     *                         post order position of the nodes. If this array
     *                         has multiple entries (i.e. multiple nodes are
     *                         selected), the node selection menu will be
     *                         positioned at the first node in this array.
     */
    SelectedNodeMenu.prototype.setSelectedNodes = function (nodeKeys) {
        // If nodeKeys includes multiple nodes, verify that all of these nodes
        // share the same name. If this _isn't_ the case, something is wrong.
        if (nodeKeys.length > 1) {
            var name = this.empress.getNodeInfo(nodeKeys[0], "name");
            for (var i = 1; i < nodeKeys.length; i++) {
                if (this.empress.getNodeInfo(nodeKeys[i], "name") !== name) {
                    throw new Error(
                        "setSelectedNodes(): keys do not represent the same " +
                            "node name."
                    );
                }
            }
        }
        /* Highlight the nodes in nodeKeys on the canvas.
         * The buffer that holds this information is formatted as
         * [x,y,r,g,b,...] where x,y are the coords of the highlighted nodes.
         */
        var highlightedNodes = [];
        for (i = 0; i < nodeKeys.length; i++) {
            var node = nodeKeys[i];
            var x = this.empress.getX(node);
            var y = this.empress.getY(node);
            highlightedNodes.push(...[x, y, 0, 1, 0]);
        }

        // send the buffer array to the drawer
        this.drawer.loadSelectedNodeBuff(highlightedNodes);

        // save the node keys to the selected node menu state machine
        this.nodeKeys = nodeKeys;
    };

    /**
     * Set the coordinates of the node menu box at the first node in nodeKeys.
     * This means that, if only a single node is selected, the menu will be
     * placed at this node's position; if multiple nodes are selected, the menu
     * will be placed at the first node's position.
     */
    SelectedNodeMenu.prototype.updateMenuPosition = function () {
        if (this.nodeKeys === null) {
            return;
        }

        var nodeToPositionAt = this.nodeKeys[0];
        // get table coords
        var x = this.empress.getX(nodeToPositionAt);
        var y = this.empress.getY(nodeToPositionAt);
        var tableLoc = this.drawer.toScreenSpace(x, y);

        // set table location. add slight offset to location so menu appears
        // next to the node instead of on top of it.
        this.box.style.left = Math.floor(tableLoc.x + 23) + "px";
        this.box.style.top = Math.floor(tableLoc.y - 43) + "px";
    };

    return SelectedNodeMenu;
});
