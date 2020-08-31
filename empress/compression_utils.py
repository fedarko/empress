# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, empress development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from collections import Counter
import pandas as pd


def remove_empty_samples_and_features(table, sample_metadata, ordination=None):
    """Removes empty samples and features from the table and sample metadata.

    This should be called *after* matching the table with the sample metadata
    and other input artifacts: we assume that the columns of the table
    DataFrame are equivalent to the indices of the sample metadata DataFrame.

    Parameters
    ----------
    table: biom.Table
        Representation of a feature table.
    sample_metadata: pd.DataFrame
        Sample metadata. The index should describe sample IDs, and the columns
        should describe sample metadata fields (e.g. "body site").
    ordination: skbio.OrdinationResults, optional
        Ordination information to show in Emperor alongside Empress. If this is
        passed, this function will check to see if any of the empty samples
        or features to be removed from the table are included in the
        ordination; if so, this will raise an error (because these empty
        items shouldn't be in the ordination in the first place).

    Returns
    -------
    filtered_table: biom.Table
        Copy of the input feature table with empty samples and features
        removed.
    filtered_sample_metadata: pd.DataFrame
        Copy of the input sample metadata with empty samples removed.

    Raises
    ------
    ValueError
        - If the input table is completely empty (i.e. all zeroes).
        - If ordination is not None, and the ordination contains empty samples
          or features.

    References
    ----------
        - Adapted from qurro._df_utils.remove_empty_samples_and_features().
    """
    orig_tbl_samples = set(table.ids())
    orig_tbl_features = set(table.ids(axis='observation'))

    # this code is equivalent to the PR below, we should update once that gets
    # merged and a newer BIOM release is publicly available
    # https://github.com/biocore/biom-format/pull/847
    filtered_table = table.copy()
    for ax in {'observation', 'sample'}:
        filtered_table = filtered_table.filter(
            table.ids(axis=ax)[table.sum(axis=ax) > 0], axis=ax,
            inplace=False)
    if filtered_table.is_empty():
        raise ValueError("All samples / features in matched table are empty.")

    # Let user know about which samples/features may have been dropped, if any.
    # Also, if we dropped any empty samples, update the sample metadata.
    filtered_sample_metadata = sample_metadata

    sample_diff = orig_tbl_samples - set(filtered_table.ids())
    if sample_diff:
        if ordination is not None:
            empty_samples_in_ord = sample_diff & set(ordination.samples.index)
            if empty_samples_in_ord:
                raise ValueError(
                    (
                        "The ordination contains samples that are empty (i.e. "
                        "all 0s) in the table. Problematic sample IDs: {}"
                    ).format(", ".join(sorted(empty_samples_in_ord)))
                )
        filtered_sample_metadata = filtered_sample_metadata.loc[
            filtered_table.ids()
        ]
        print("Removed {} empty sample(s).".format(len(sample_diff)))

    feature_diff = orig_tbl_features - \
        set(filtered_table.ids(axis='observation'))
    if feature_diff:
        if ordination is not None and ordination.features is not None:
            empty_feats_in_ord = feature_diff & set(ordination.features.index)
            if empty_feats_in_ord:
                raise ValueError(
                    (
                        "The ordination contains features that are empty "
                        "(i.e. all 0s) in the table. Problematic feature IDs: "
                        "{}"
                    ).format(", ".join(sorted(empty_feats_in_ord)))
                )
        print("Removed {} empty feature(s).".format(len(feature_diff)))

    return filtered_table, filtered_sample_metadata


def compress_table(table):
    """Converts a feature table to a space-saving format.

    Parameters
    ----------
    table: biom.Table
        Representation of a feature table.  It is assumed that empty samples /
        features have already been removed from the table.

    Returns
    -------
    (s_ids, f_ids, s_ids_to_indices, f_ids_to_indices, compressed_table)
        s_ids: list
            List of the sample IDs in the table.
        f_ids: list
            List of the feature IDs in the table, analogous to s_ids.
        s_ids_to_indices: dict
            Inverse of s_ids: this maps sample IDs to their indices in s_ids.
            "Indices" refers to a feature or sample's 0-based position in f_ids
            or s_ids, respectively.
        f_ids_to_indices: dict
            Inverse of f_ids: this maps feature IDs to their indices in f_ids,
            analogous to s_ids_to_indices.
        compressed_table: list
            Two-dimensional list. The "outer list" is of length len(s_ids).
            Each position i within this outer list holds an "inner list" of
            arbitrary (but within the range [1, len(f_ids)]) length.
            The i-th inner list contains the feature indices of the
            features present (i.e. at any abundance > 0) within the
            sample with index i. Each inner list is sorted in ascending order.

    References
    ----------
        - Inspired by redbiom and Qurro's JSON data models.
    """
    feature_ids = table.ids(axis='observation')
    sample_ids = table.ids()

    f_ids_to_indices = {fid: idx for idx, fid in enumerate(feature_ids)}
    s_ids_to_indices = {sid: idx for idx, sid in enumerate(sample_ids)}

    compressed_table = []
    for vec in table.iter_data(axis='sample', dense=False):
        compressed_table.append([int(i) for i in vec.indices])

    return (
        list(sample_ids), list(feature_ids), s_ids_to_indices,
        f_ids_to_indices, compressed_table
    )


def compress_recurring_md_vals(str_metadata_df, output_type="list"):
    """Replaces recurring values in a metadata DF with integers.

    Parameters
    ----------
    str_metadata_df: pd.DataFrame
        Sample or feature metadata. Index should describe samples (or
        features), and columns should be metadata fields. All of the values in
        this DataFrame must be stored as strings (even numbers should be stored
        as strings), so that compressed values can be interpreted correctly.

    Returns
    -------
    (recurring_vals, compressed_md)
        recurring_vals: list
            List of "recurring values" that are used at least twice in the
            metadata, sorted in descending order of frequency in the metadata
            (with ties broken arbitrarily).
        compressed_md: list
            Two-dimensional list. The "outer list" is of length
            len(str_metadata_df.index). Each position i within this outer
            list holds an "inner list" of length len(str_metadata_df.columns).
            The c-th value of the i-th inner list contains either:
             -An integer, in which case the value referred to is located at
              this integer's position in recurring_vals (using 0-indexing)
             -A string value
            In either case, the value referred to is the metadata value in
            column c for the sample with index i.
    """
    num_cols = len(str_metadata_df.columns)
    orig_idx = list(str_metadata_df.index)

    # Generate a numpy ndarray of all metadata values
    vals_flattened = str_metadata_df.values.flatten()

    # Count how many times each metadata value occurs, and (with most_common())
    # sort these values in descending order of frequency. Note that the "sort
    # in descending order" step isn't technically needed; although having the
    # most frequent value be at position 0 in recurring_vals, etc. makes
    # interpreting the compression results easier (and might save some space
    # for large-enough datasets since writing "0" 1000 times takes up less
    # characters than writing "1000" 1000 times), it's probably not a big deal.
    vals_and_counts = Counter(vals_flattened).most_common()

    # We temporarily store a mapping of each recurring value to its position in
    # recurring_vals, to make populating compressed_md easier.
    val2cvidx = {}
    recurring_vals = []
    for (val, ct) in vals_and_counts:
        if ct > 1:
            val2cvidx[val] = len(recurring_vals)
            recurring_vals.append(val)
        else:
            # Since vals_and_counts is in descending order, once we get to a
            # value that only occurs once (i.e. this value is not "recurring")
            # we can stop iterating through vals_and_counts.
            break

    if output_type == "list":
        # compressed_md will be a 2-D list representation of the input metadata
        # with each recurring value replaced with an integer pointing to a
        # position in recurring_vals.
        compressed_md = []
        c = 0
        for val in vals_flattened:
            if c % num_cols == 0:
                compressed_md.append([])
            # Add this value (either a number in recurring_vals or the actual
            # string value) to the last list in compressed_md, corresponding to
            # the current sample or feature
            if val in val2cvidx:
                compressed_md[-1].append(val2cvidx[val])
            else:
                compressed_md[-1].append(val)
            c += 1
    elif output_type == "dict":
        # compressed_md will be a dict representation of the input metadata,
        # where each value is a list of length num_cols (with recurring values
        # replaced as above).
        compressed_md = {}
        c = 0
        r = 0
        curr_feature_name = None
        for val in vals_flattened:
            if c % num_cols == 0:
                compressed_md[orig_idx[r]] = []
                curr_feature_name = orig_idx[r]
                r += 1
            if val in val2cvidx:
                compressed_md[curr_feature_name].append(val2cvidx[val])
            else:
                compressed_md[curr_feature_name].append(val)
            c += 1
    else:
        raise ValueError('Unrecognized output_type; must be "list" or "dict"')

    return recurring_vals, compressed_md


def compress_sample_metadata(s_ids_to_indices, metadata):
    """Converts a sample metadata DataFrame to a space-saving format.

    Parameters
    ----------
    s_ids_to_indices: dict
        Maps sample IDs (strings) to 0-based indices in an existing list of
        sample IDs. In practice, this should just be the "s_ids_to_indices"
        output from compress_table().
    metadata: pd.DataFrame
        Sample metadata. The index should describe sample IDs, and the columns
        should describe sample metadata fields (e.g. "body site").
        The sample IDs in the index should match one-to-one with the keys in
        s_ids_to_indices.

    Returns
    -------
    (sm_cols, recurring_vals, compressed_sm)
        sm_cols: list
            List of the sample metadata column names, all converted to strings.
        recurring_vals: list
            List of "recurring values" that are used at least twice in the
            sample metadata, sorted in descending order of frequency in the
            sample metadata (with ties broken arbitrarily).
        compressed_sm: list
            Two-dimensional list. The "outer list" is of length
            len(s_ids_to_indices.keys()). Each position i within this outer
            list holds an "inner list" of length len(sm_cols).
            The c-th value of the i-th inner list contains either:
             -An integer, in which case the value referred to is located at
              this integer's position in recurring_vals (using 0-indexing)
             -A string value
            In either case, the value referred to is the sample metadata value
            in column c for the sample with index i.

    Raises
    ------
    ValueError
        - If the metadata's index and the keys of s_ids_to_indices do not
          contain the exact same elements.
        - If the values of s_ids_to_indices are invalid: that is, if sorting
          the values in ascending order does not produce a list of
          [0, 1, 2, 3, ..., len(s_ids_to_indices.keys())].

    References
    ----------
        - Inspired by redbiom and Qurro's JSON data models.
    """
    sample_ids = s_ids_to_indices.keys()
    # NOTE: I think that identically-named samples or metadata columns will
    # break this check, but I also think that we can assume by this point that
    # the data is at least that sane. (Checking that should be a responsibility
    # for earlier in the program.)
    if set(sample_ids) != set(metadata.index):
        raise ValueError(
            "The sample IDs in the metadata's index and s_ids_to_indices are "
            "not identical."
        )

    if sorted(s_ids_to_indices.values()) != list(range(len(sample_ids))):
        raise ValueError("Indices (values) of s_ids_to_indices are invalid.")

    # Rename sample IDs to indices in the metadata
    indexed_metadata = metadata.rename(index=s_ids_to_indices)

    # Sort the metadata's rows by the sample indices
    sorted_i_metadata = indexed_metadata.sort_index(
        axis="index", ascending=True
    )

    # Convert all of the metadata values to strings
    str_s_i_metadata = sorted_i_metadata.astype(str)

    # Also convert the metadata columns to strings
    sm_cols = [str(c) for c in str_s_i_metadata.columns]

    # Compress recurring values and convert the metadata to a 2-D list
    recurring_vals, compressed_sm = compress_recurring_md_vals(
        str_s_i_metadata
    )

    return sm_cols, recurring_vals, compressed_sm


def compress_feature_metadata(tip_metadata, int_metadata, name2treepos):
    """Converts tip/internal node metadata DataFrames to a space-saving format.

    Note that the columns of tip_metadata and int_metadata should be identical,
    even if the feature metadata only describes tip or internal nodes. (In that
    case, then the other feature metadata parameter should still be a DataFrame
    -- albeit an empty one, with no feature names in its index.) The only case
    in which either of the _metadata parameters should be None is if there was
    no feature metadata at all.

    Parameters
    ----------
    tip_metadata: pd.DataFrame or None
        Metadata for tip nodes. If not None, the index should describe node
        names, and the columns should describe feature metadata fields.
    int_metadata: pd.DataFrame or None
        Metadata for internal nodes. If not None, the index should describe
        node names, and the columns should describe feature metadata fields.
    name2treepos: dict
        Maps node names to a list of corresponding postorder position(s) in the
        tree. (Tip names should only map to one position, while internal node
        names can map to an arbitrarily large number of positions.)

    Returns
    -------
    (fm_cols, recurring_vals, compressed_tm, compressed_im)
        fm_cols: list
            List of the feature metadata column names, all converted to
            strings. If both input DFs are None, this will be {}.
        recurring_vals: list
            List of "recurring values" that are used at least twice in either
            the tip or internal node metadata, sorted in descending order of
            frequency in the metadata (with ties broken arbitrarily). (Note
            that these are computed by looking at both metadata DFs
            at once: so if a given value is used once in the tip metadata and
            once in the internal node metadata then it'll still get included in
            this list.
        compressed_tm: list
            Two-dimensional list representation of the tip metadata. Along with
            recurring values being replaced with their position in
            recurring_vals, the DF's indices (which were previously tip node
            names in the tree) will be replaced with each tip's postorder
            position (based on name2treepos). If tip_metadata was empty, or if
            both input DFs were None, this will be {}.
        compressed_im: list
            Two-dimensional list representation of the internal node metadata.
            Recurring values are replaced as described above for compressed_tm,
            and the DF's indices will similarly be replaced with each internal
            node's postorder position based on name2treepos. Rows where the
            index (a.k.a. internal node name) maps to multiple postorder
            positions in name2treepos will be duplicated, such that each tree
            position with internal node metadata will have its own row. (This
            could be made more efficient, but it simplifies things a lot.)
            If tip_metadata was empty, or if both input DFs were None, this
            will be {}.

    Raises
    ------
    ValueError
        - If only one of tip_metadata and int_metadata is None.
        - If the columns of tip_metadata are not identical to the columns of
          int_metadata.
        - If both the tip and internal node metadata DataFrames are empty.
        - If any of the node names in tip_metadata maps to multiple node
          positions in name2treepos.

    KeyError
        - If any of the node names (indices) in tip_metadata or int_metadata
          are not present as keys in name2treepos.

    References
    ----------
        - Inspired by redbiom and Qurro's JSON data models.
    """
    # If the user didn't pass in any feature metadata, we'll get to this block
    if tip_metadata is None and int_metadata is None:
        return [], [], {}, {}

    # *This* should never happen. If it did, it's a sign that this function is
    # being misused. (The ^ is a logical XOR; see
    # https://stackoverflow.com/a/432844/10730311.)
    if (tip_metadata is None) ^ (int_metadata is None):
        raise ValueError(
            "Only one of tip & int. node feature metadata is None."
        )

    # Verify that columns match up btwn. tip and internal node metadata
    if not tip_metadata.columns.equals(int_metadata.columns):
        raise ValueError("Tip & int. node feature metadata columns differ.")

    # Verify that at least one feature metadata entry exists (since at this
    # point we know that there should be at least *some* feature metadata)
    if tip_metadata.empty and int_metadata.empty:
        raise ValueError("Both tip & int. node feature metadata are empty.")

    # Make sure that all columns and values are passed to the JS as strings
    fm_cols = [str(c) for c in tip_metadata.columns]
    str_tm = tip_metadata.astype(str)
    str_im = int_metadata.astype(str)

    # Combine the tip / internal node metadata so that they can share recurring
    # values
    combined_str_fm = pd.concat([str_tm, str_im])

    # Compress recurring values
    recurring_vals, compressed_fm = compress_recurring_md_vals(
        combined_str_fm, output_type="dict"
    )

    # Split up the compressed feature metadata back into tip and internal node
    # metadata. Also, while we're going through the metadata, replace node
    # names with the corresponding indices in name2treepos, and duplicate
    # internal node metadata.
    compressed_tm = {}
    compressed_im = {}
    for name in compressed_fm:
        try:
            treepositions = name2treepos[name]
        except KeyError:
            raise KeyError(
                "Node name {} is not present in name2treepos.".format(name)
            )
        if name in tip_metadata.index:
            if len(treepositions) > 1:
                raise ValueError("Tip node name is shared by multiple nodes.")
            compressed_tm[treepositions[0]] = compressed_fm[name]
        else:
            for pos in treepositions:
                compressed_im[pos] = compressed_fm[name]

    return fm_cols, recurring_vals, compressed_tm, compressed_im
