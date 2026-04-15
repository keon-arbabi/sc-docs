.. _typedefs:

Type aliases
============

.. currentmodule:: single_cell

SingleCellColumn
----------------

``str | pl.Expr | pl.Series | np.ndarray | Callable[[SingleCell], pl.Series | np.ndarray]``

A flexible column specification. Can be a column name (string), a Polars
expression, a Polars Series, a 1D NumPy array, or a function that takes a
SingleCell dataset and returns a Series or array.

Used by :meth:`SingleCell.qc`, :meth:`SingleCell.find_doublets`,
:meth:`SingleCell.pseudobulk`, :meth:`SingleCell.hvg`,
:meth:`SingleCell.normalize`, and others.

UnsDict
-------

``dict[str, str | int | float | bool | np.ndarray | 'UnsDict']``

A nested dictionary of unstructured metadata. Values can be scalars
(strings, numbers, Booleans), NumPy arrays, or nested dictionaries
of the same types. Stored in :attr:`SingleCell.uns`.
