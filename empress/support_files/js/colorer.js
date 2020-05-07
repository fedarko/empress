define(["chroma", "underscore", "util"], function (chroma, _, util) {
    /**
     * @class Colorer
     *
     * Creates a color object that will map values to colors from a pre-defined
     * color map.
     *
     * @param{Object} color The color map to draw colors from.
     *                      This should be an id in Colorer.__Colormaps.
     * @param{Array} values The values in a metadata field for which colors
     *                      will be generated.
     *
     * @return{Colorer}
     * constructs Colorer
     */
    function Colorer(color, values) {
        // Remove duplicate values and sort the values sanely
        var sortedValues = util.naturalSort(_.uniq(values));

        this.__valueToColor = {};

        // Figure out what "type" of color map has been requested (discrete vs.
        // sequential / diverging). Will inform how we assign values to colors.
        var type;
        var name;
        var colormap;
        for (var c = 0; c < Colorer.__Colormaps.length; c++) {
            colormap = Colorer.__Colormaps[c];
            if (colormap.id === color) {
                type = colormap.type;
                name = colormap.name;
                break;
            }
        }
        if (type === Colorer.DISCRETE) {
            var palette;
            if (color === Colorer.__QIIME_COLOR) {
                palette = Colorer.__qiimeDiscrete;
            } else {
                palette = chroma.brewer[color];
            }
            for (var i = 0; i < sortedValues.length; i++) {
                var modIndex = i % palette.length;
                this.__valueToColor[sortedValues[i]] = palette[modIndex];
            }
        } else if (type === Colorer.SEQUENTIAL || type === Colorer.DIVERGING) {
            var map = chroma.brewer[color];

            // set up a closure so we can update this.__valueToColor within the
            // functions below
            var thisColorer = this;

            // Get list of only numeric values, and see if it's possible for us
            // to do scaling with these values or not
            var split = util.splitNumericValues(sortedValues);

            if (split.numeric.length < 2) {
                // We don't have enough numeric values to do any sort of
                // "scaling," so instead we just color everything gray
                _.each(values, function (val) {
                    thisColorer.__valueToColor[val] = Colorer.NANCOLOR;
                });
                alert(
                    "The selected color map (" +
                        name +
                        ") cannot be used " +
                        "with this field. Continuous coloration requires at " +
                        "least 2 numeric values in the field."
                );
            } else {
                // We can do scaling! Nice.
                // convert objects to numbers so we can map them to a color
                numbers = _.map(split.numeric, parseFloat);
                min = _.min(numbers);
                max = _.max(numbers);

                var interpolator = chroma.scale(map).domain([min, max]);

                // Color all the numeric values
                _.each(split.numeric, function (element) {
                    thisColorer.__valueToColor[element] = interpolator(
                        +element
                    );
                });
                // Gray out non-numeric values
                _.each(split.nonNumeric, function (element) {
                    thisColorer.__valueToColor[element] = Colorer.NANCOLOR;
                });
            }
        } else {
            throw new Error("Invalid color " + color + " specified");
        }
    }

    /**
     * Returns an rgb array with values in the range of [0, 1].
     *
     * @param{String} value A value that was present in the values used to
     *                      construct this Colorer object.
     *
     * @return{Object} An rgb array
     */
    Colorer.prototype.getColorRGB = function (value) {
        // the slice() strips off the opacity element, which causes problems
        // with Empress' drawing code
        return chroma(this.__valueToColor[value]).gl().slice(0, 3);
    };

    /**
     * Returns an rgb hex string.
     *
     * @param{String} value A value that was present in the values used to
     *                      construct this Colorer object.
     *
     * @return{Object} An rgb hex string
     */
    Colorer.prototype.getColorHex = function (value) {
        return this.__valueToColor[value];
    };

    /**
     * Adds all available color maps to the select object.
     *
     * @param{Object} sel The select object to add color map options to.
     * @classmethod
     */
    Colorer.addColorsToSelect = function (sel) {
        // The color map selector
        for (var i = 0; i < Colorer.__Colormaps.length; i++) {
            var map = Colorer.__Colormaps[i];
            var opt = document.createElement("option");
            opt.innerHTML = map.name;
            opt.value = map.id;

            if (map.type == "Header") {
                opt.disabled = true;
            }
            sel.appendChild(opt);
        }
    };

    Colorer.DISCRETE = "Discrete";
    Colorer.SEQUENTIAL = "Sequential";
    Colorer.DIVERGING = "Diverging";
    Colorer.HEADER = "Header";

    // taken from the qiime/colors.py module; a total of 24 colors
    /** @private */
    Colorer.__QIIME_COLOR = "discrete-coloring-qiime";
    Colorer.__qiimeDiscrete = [
        "#ff0000",
        "#0000ff",
        "#f27304",
        "#008000",
        "#91278d",
        "#ffff00",
        "#7cecf4",
        "#f49ac2",
        "#5da09e",
        "#6b440b",
        "#808080",
        "#f79679",
        "#7da9d8",
        "#fcc688",
        "#80c99b",
        "#a287bf",
        "#fff899",
        "#c49c6b",
        "#c0c0c0",
        "#ed008a",
        "#00b6ff",
        "#a54700",
        "#808000",
        "#008080",
    ];

    // This is also the default "nanColor" for Emperor. (We could make this
    // configurable if desired.)
    Colorer.NANCOLOR = "#64655d";

    // Used to create color select option and chroma.brewer
    //Modified from:
    //https://github.com/biocore/emperor/blob/
    //     027aa16f1dcf9536cd2dd9c9800ece5fc359ecbc/emperor/
    //     support_files/js/color-view-controller.js#L573-L613
    Colorer.__Colormaps = [
        { name: "-- Discrete --", type: Colorer.HEADER },
        {
            id: "discrete-coloring-qiime",
            name: "Classic QIIME Colors",
            type: Colorer.DISCRETE,
        },
        { id: "Paired", name: "Paired", type: Colorer.DISCRETE },
        { id: "Accent", name: "Accent", type: Colorer.DISCRETE },
        { id: "Dark2", name: "Dark", type: Colorer.DISCRETE },
        { id: "Set1", name: "Set1", type: Colorer.DISCRETE },
        { id: "Set2", name: "Set2", type: Colorer.DISCRETE },
        { id: "Set3", name: "Set3", type: Colorer.DISCRETE },
        { id: "Pastel1", name: "Pastel1", type: Colorer.DISCRETE },
        { id: "Pastel2", name: "Pastel2", type: Colorer.DISCRETE },

        { name: "-- Sequential --", type: Colorer.HEADER },
        { id: "Viridis", name: "Viridis", type: Colorer.SEQUENTIAL },
        { id: "Reds", name: "Reds", type: Colorer.SEQUENTIAL },
        { id: "RdPu", name: "Red-Purple", type: Colorer.SEQUENTIAL },
        { id: "Oranges", name: "Oranges", type: Colorer.SEQUENTIAL },
        { id: "OrRd", name: "Orange-Red", type: Colorer.SEQUENTIAL },
        { id: "YlOrBr", name: "Yellow-Orange-Brown", type: Colorer.SEQUENTIAL },
        { id: "YlOrRd", name: "Yellow-Orange-Red", type: Colorer.SEQUENTIAL },
        { id: "YlGn", name: "Yellow-Green", type: Colorer.SEQUENTIAL },
        { id: "YlGnBu", name: "Yellow-Green-Blue", type: Colorer.SEQUENTIAL },
        { id: "Greens", name: "Greens", type: Colorer.SEQUENTIAL },
        { id: "GnBu", name: "Green-Blue", type: Colorer.SEQUENTIAL },
        { id: "Blues", name: "Blues", type: Colorer.SEQUENTIAL },
        { id: "BuGn", name: "Blue-Green", type: Colorer.SEQUENTIAL },
        { id: "BuPu", name: "Blue-Purple", type: Colorer.SEQUENTIAL },
        { id: "Purples", name: "Purples", type: Colorer.SEQUENTIAL },
        { id: "PuRd", name: "Purple-Red", type: Colorer.SEQUENTIAL },
        { id: "PuBuGn", name: "Purple-Blue-Green", type: Colorer.SEQUENTIAL },
        { id: "Greys", name: "Greys", type: Colorer.SEQUENTIAL },

        { name: "-- Diverging --", type: Colorer.HEADER },
        { id: "Spectral", name: "Spectral", type: Colorer.DIVERGING },
        { id: "RdBu", name: "Red-Blue", type: Colorer.DIVERGING },
        { id: "RdYlGn", name: "Red-Yellow-Green", type: Colorer.DIVERGING },
        { id: "RdYlBu", name: "Red-Yellow-Blue", type: Colorer.DIVERGING },
        { id: "RdGy", name: "Red-Grey", type: Colorer.DIVERGING },
        { id: "PiYG", name: "Pink-Yellow-Green", type: Colorer.DIVERGING },
        { id: "BrBG", name: "Brown-Blue-Green", type: Colorer.DIVERGING },
        { id: "PuOr", name: "Purple-Orange", type: Colorer.DIVERGING },
        { id: "PRGn", name: "Purple-Green", type: Colorer.DIVERGING },
    ];
    return Colorer;
});
