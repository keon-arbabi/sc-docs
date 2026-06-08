:orphan:

.. _typedefs:

Type aliases
============

.. currentmodule:: brisc

SingleCellColumn
----------------

``str | pl.Expr | pl.Series | np.ndarray | Callable[[SingleCell], pl.Series | np.ndarray]``

A flexible column specification. Can be a column name (string), a Polars
expression, a Polars Series, a 1D NumPy array, or a function that takes a
SingleCell dataset and returns a Series or array.

Used by :meth:`SingleCell.qc`, :meth:`SingleCell.find_doublets`,
:meth:`SingleCell.pseudobulk`, :meth:`SingleCell.hvg`,
:meth:`SingleCell.normalize`, and others.

Scalar
------

``str | int | float | Decimal | date | time | datetime | timedelta | bool | bytes | Expr | Series | Iterable[Scalar]``

Any primitive literal or sequence thereof. Accepted wherever a single
value (or a 1-D sequence of values) is expected -- filter thresholds,
default fills, single-column inputs, etc.

UnsDict
-------

``dict[str, str | int | float | bool | np.ndarray | 'UnsDict']``

A nested dictionary of unstructured metadata. Values can be scalars
(strings, numbers, Booleans), NumPy arrays, or nested dictionaries
of the same types. Stored in :attr:`SingleCell.uns`.
