from skbio import TreeNode
import numpy as np


class Tree(TreeNode):
    """
    Attributes
    ----------
    length
    leafcount
    height
    depth

    Notes
    -----
    `length` refers to the branch length of a node to its parent.
    `leafcount` is the number of tips within a subtree. `height` refers
    to the longest path from root to the deepst leaf in that subtree.
    `depth` is the number of nodes found in the longest path.
    """

    def __init__(self, use_lengths=False, **kwargs):
        """ Constructs a Dendrogram object for visualization.

        Parameters
        ----------
        use_lengths: bool
            Specifies if the branch lengths should be included in the
            resulting visualization (default True).

        Returns
        -------

        """
        super().__init__(**kwargs)
        self.childRem = -1

    @classmethod
    def from_tree(cls, tree, use_lengths=True):
        """ Creates an Tree object from a skbio tree.

        Parameters
        ----------
        tree : skbio.TreeNode
            Input skbio tree
        use_lengths: Boolean
            Specify if the branch length should be incorporated into
            the geometry calculations for visualization.
        Returns
        -------
        Tree

        """
        if tree.count() <= 1:
            raise ValueError("Tree must contain at least 2 nodes.")

        for n in tree.postorder():
            n.__class__ = Tree
            n.tip_count = 0

        tree.update_geometry(use_lengths)
        return tree

    def update_geometry(self, use_lengths, depth=None):
        """Calculate tree node attributes such as height and depth.

        Parameters
        ----------
        use_lengths: bool
            Specify if the branch length should be incorporated into
            the geometry calculations for visualization.
        depth: int
            The number of nodes in the longest path from root to leaf.
            This is agnostic to scale and orientation.

        """
        # i = 0
        for node in self.postorder():
            # node.name = i
            # i += 1
            if node.length is None or not use_lengths:
                if not use_lengths:
                    if node.is_tip():
                        node.length = 5
                    else:
                        node.length = 1
                else:
                    node.length = 0

            node.depth = (depth or 0) + node.length

            children = node.children
            if children:
                node.height = max([c.height for c in children]) + node.length
                node.leafcount = sum([c.leafcount for c in children])

            else:
                node.height = node.length
                node.leafcount = 1

    def coords(self, height, width):
        """ Computes the coordinates of nodes to be rendered in plot.

        This runs multiple layout algorithms and saves all of the resulting
        coordinates for each node, so that layout algorithms can be rapidly
        toggled between in the JS interface.

        Also adds on .highestchildyr and .lowestchildyr attributes to internal
        nodes so that vertical bars for these nodes can be drawn in the
        rectangular layout.

        Parameters
        ----------
        height : int
            The height of the canvas.
        width : int
            The width of the canvas.

        Returns
        -------
        layout_to_coordsuffix : dict
        default_layout : str
        """

        layout_to_coordsuffix = {}
        layout_algs = (
            self.layout_unrooted,
            self.layout_rectangular,
        )
        # We set the default layout to whatever the first layout in
        # layout_algs is, but this behavior is of course modifiable
        default_layout = None
        for alg in layout_algs:
            name, suffix = alg(width, height)
            layout_to_coordsuffix[name] = suffix
            self.alter_coordinates_relative_to_root(suffix)
            if default_layout is None:
                default_layout = name

        # Determine highest and lowest child y-position for internal nodes in
        # the rectangular layout; used to draw vertical lines for these nodes.
        #
        # NOTE / TODO: This will have the effect of drawing vertical lines even
        # for nodes with only 1 child -- in this case lowestchildyr ==
        # highestchildyr for this node, so all of the stuff drawn in WebGL for
        # this vertical line shouldn't show up. I don't think this should cause
        # any problems, but it may be worth detecting these cases and not
        # drawing vertical lines for them in the future.
        for n in self.preorder():
            if not n.is_tip():
                child_y_coords = [c.yr for c in n.children]
                n.highestchildyr = max(child_y_coords)
                n.lowestchildyr = min(child_y_coords)

        return layout_to_coordsuffix, default_layout

    def alter_coordinates_relative_to_root(self, suffix):
        """ Subtracts the root node's x- and y- coords from all nodes' coords.

        This was previously done within coords(), but I moved it here so that
        this logic can be used after arbitrary layout computations.

        Parameters
        ----------
        suffix : str
            The suffix of the x- and y-coordinates to adjust.

            For example, this is "2" for the unrooted layout since coordinates
            are stored in the x2 and y2 attributes for every node; and it's "r"
            for the rectangular layout since the coordinate attributes are now
            xr and yr.
        """

        xname = "x" + suffix
        yname = "y" + suffix

        centerX = getattr(self, xname)
        centerY = getattr(self, yname)

        for node in self.postorder():
            # This code might look sort of intimidating, but it's really just
            # another way to write out:
            #     node.x2 = node.x2 - centerX
            #     node.y2 = node.y2 - centerY
            # ...when we don't know what "x2" or "y2" will be named beforehand.
            setattr(node, xname, getattr(node, xname) - centerX)
            setattr(node, yname, getattr(node, yname) - centerY)

    def layout_circular(self, width, height):
        """ Circular layout version of the rectangular layout.

        Currently incomplete, but in theory this should work analogously to the
        rectangular layout: start at the root (with available angle ranges [0,
        2pi], and go down through the tree's nodes in a preorder(?) traversal,
        dividing the angles as needed based on the number of leaves of a child
        node.

        Parameters
        ----------
        width : float
            width of the canvas
        height : float
            height of the canvas

        References
        ----------
        https://github.com/qiime/Topiary-Explorer/blob/master/src/topiaryexplorer/TreeVis.java
            Description above + the implementation of this algorithm
            derived from the Polar/Radial layout algorithm code.
        """
        pass

    def layout_rectangular(self, width, height):
        """ Rectangular layout.

        In this sort of layout, each tip has a distinct y-position, and parent
        y-positions are centered over their descendant tips' positions.
        x-positions are computed based on nodes' branch lengths.

        Following this algorithm, nodes' unrooted layout coordinates are
        accessible at [node].xr and [node].yr.

        For a simple tree, this layout should look something like:
                 __
             ___|
         ___|   |__
        |   |___
        |    ___
        |___|
            |___

        Parameters
        ----------
        width : float
            width of the canvas
        height : float
            height of the canvas

        References
        ----------
        https://rachel53461.wordpress.com/2014/04/20/algorithm-for-drawing-trees/
            Clear explanation of Reingold-Tilford that I used a lot
        https://github.com/qiime/Topiary-Explorer/blob/master/src/topiaryexplorer/TreeVis.java
            Derived from the "Rectangular" layout algorithm code.
        """
        # NOTE: This doesn't draw a horizontal line leading to the root "node"
        # of the graph. See https://github.com/biocore/empress/issues/141 for
        # context.
        max_width = 0
        max_height = 0
        prev_y = 0
        for n in self.postorder():
            if n.is_tip():
                n.yr = prev_y
                prev_y += 1
                if n.yr > max_height:
                    max_height = n.yr
            else:
                # Center internal nodes above their children
                # We could also center them above their tips, but (IMO) this
                # looks better ;)
                n.yr = sum([c.yr for c in n.children]) / len(n.children)

        self.xr = 0
        for n in self.preorder(include_self=False):
            n.xr = n.parent.xr + n.length
            if n.xr > max_width:
                max_width = n.xr

        if max_width == 0:
            # The only case in which this should happen is if all nodes in the
            # tree have a branch length of 0.
            #
            # We could theoretically allow this case (especially if/when we add
            # in functionality to ignore branch lengths and just draw every
            # node with length of "1" or something), but to me this seems like
            # the sort of thing we should flag since that isn't currently
            # implemented.
            raise ValueError(
                "All nodes in the tree have branch lengths of 0."
            )
        x_scaling_factor = width / max_width

        if max_height > 0:
            # Having a max_height of 0 could actually happen, in the funky case
            # where the entire tree is a straight line (e.g. A -> B -> C). In
            # this case our "rectangular layout" drawing places all nodes on
            # the same y-coordinate (0), resulting in max_height = 0.
            # ... So, that's why we only do y-scaling if this *isn't* the case.
            y_scaling_factor = height / max_height
        else:
            # Since this will be multiplied by 0 for every node, we can set
            # this to any real number and get the intended "effect" of keeping
            # every node's y-coordinate at 0.
            y_scaling_factor = 1

        for n in self.preorder():
            n.xr *= x_scaling_factor
            n.yr *= y_scaling_factor

        # Now we have the layout! In the JS we'll need to draw each internal
        # node as a vertical line ranging from its lowest child y-position to
        # its highest child y-position, and then draw horizontal lines from
        # this line to all of its child nodes (where the length of the
        # horizontal line is proportional to the node length in question).
        return "Rectangular", "r"

    def layout_unrooted(self, width, height):
        """ Find best scaling factor for fitting the tree in the figure.
        This method will find the best orientation and scaling possible to
        fit the tree within the dimensions specified by width and height, using
        an unrooted layout algorithm.

        Following this algorithm, nodes' unrooted layout coordinates are
        accessible at [node].x2 and [node].y2.

        Parameters
        ----------
        width : float
            width of the canvas
        height : float
            height of the canvas

        Returns
        -------
        best_scaling : float
            largest scaling factor in which the tree can fit in the canvas.

        Notes
        -----

        """
        # Recall that 360 degrees is equal to (2 * pi) radians.
        # You can think of this variable as "the maximum angle we can 'give' to
        # each leaf of the tree".
        angle = (2 * np.pi) / self.leafcount

        best_scale = 0
        for i in range(60):
            direction = i / 60.0 * np.pi

            (max_x, min_x, max_y, min_y) = self.update_unrooted_coords(
                1.0, 0, 0, direction, angle)

            x_diff = max_x - min_x
            width_min = 0
            if x_diff != 0:
                width_min = float(width) / x_diff
            y_diff = max_y - min_y
            height_min = 0
            if y_diff != 0:
                height_min = float(height) / y_diff
            scale = min(width_min, height_min)
            scale *= 0.95  # extra margin for labels
            if scale >= best_scale:
                best_scale = scale
                mid_x = width / 2 - ((max_x + min_x) / 2) * scale
                mid_y = height / 2 - ((max_y + min_y) / 2) * scale
                best_args = (scale, mid_x, mid_y, direction, angle)

        self.update_unrooted_coords(*best_args)
        return "Unrooted", "2"

    def update_unrooted_coords(self, s, x1, y1, a, da):
        """ Update x, y coordinates of tree nodes in canvas.

        This function will update the x1, y1, x2, y2, and angle attributes
        for all of the nodes within the tree. Note that (once the unrooted
        layout has finished) all that is really used are the x2 and y2
        attributes.

        In a server-based version of Empress, this could be applied when
        the tree becomes modified (i.e. pruning or collapsing) and the
        resulting coordinates would be modified to reflect the changes
        to the tree structure. (In practice, we just run this once on the
        Python side of things in order to precompute the layout.)

        Parameters
        ----------
        s : float
            scaling
        x1 : float
            x midpoint
        y1 : float
            y midpoint
        a : float
            angle (degrees)
        da : float
            angle resolution (degrees)

        Returns
        -------
        points : list of tuple
            2D coordinates of all of the nodes.
        """

        max_x = float('-inf')
        min_x = float('inf')
        max_y = float('-inf')
        min_y = float('inf')

        # calculates self coords/angle
        # Constant angle algorithm.  Should add maximum daylight step.
        x2 = x1 + self.length * s * np.sin(a)
        y2 = y1 + self.length * s * np.cos(a)
        (self.x1, self.y1, self.x2, self.y2, self.angle) = (x1, y1, x2, y2,
                                                            a)
        nodes = [node for node in self.postorder(include_self=False)]
        nodes.reverse()
        # for node in self.preorder(include_self=False):
        for node in nodes:
            x1 = node.parent.x2
            y1 = node.parent.y2

            # init a
            a = node.parent.angle

            # same modify across nodes
            a = a - node.parent.leafcount * da / 2

            # check for conditional higher order
            for sib in node.parent.children:
                if sib != node:
                    a += sib.leafcount * da
                else:
                    a += (node.leafcount * da) / 2
                    break

            # Constant angle algorithm.  Should add maximum daylight step.
            x2 = x1 + node.length * s * np.sin(a)
            y2 = y1 + node.length * s * np.cos(a)
            (node.x1, node.y1, node.x2, node.y2, node.angle) = (x1, y1, x2,
                                                                y2, a)

            max_x, min_x = max(max_x, x2), min(min_x, x2)
            max_y, min_y = max(max_y, y2), min(min_y, y2)

        return (max_x, min_x, max_y, min_y)
