Pseudobulk
==========

.. toctree::
   :maxdepth: 2
   :hidden:

   constructor
   io
   properties
   data_access
   dictionary_interface
   manipulation
   transformation
   analysis
   utility

.. currentmodule:: single_cell

.. autoclass:: Pseudobulk
   :no-members:

Constructor
-----------

.. autosummary::
   :nosignatures:

   Pseudobulk.__init__

I/O
---

.. autosummary::
   :nosignatures:

   Pseudobulk.save

Properties
----------

.. autosummary::
   :nosignatures:

   Pseudobulk.X
   Pseudobulk.obs
   Pseudobulk.var
   Pseudobulk.obs_names
   Pseudobulk.var_names
   Pseudobulk.num_threads
   Pseudobulk.shape

Data access
-----------

.. autosummary::
   :nosignatures:

   Pseudobulk.sample
   Pseudobulk.gene

Dictionary interface
--------------------

.. autosummary::
   :nosignatures:

   Pseudobulk.keys
   Pseudobulk.values
   Pseudobulk.items
   Pseudobulk.iter_X
   Pseudobulk.iter_obs
   Pseudobulk.iter_var
   Pseudobulk.__contains__
   Pseudobulk.__or__
   Pseudobulk.__eq__

Manipulation
------------

.. autosummary::
   :nosignatures:

   Pseudobulk.set_obs_names
   Pseudobulk.set_var_names
   Pseudobulk.set_num_threads
   Pseudobulk.filter_obs
   Pseudobulk.filter_var
   Pseudobulk.select_obs
   Pseudobulk.select_var
   Pseudobulk.select_cell_types
   Pseudobulk.with_columns_obs
   Pseudobulk.with_columns_var
   Pseudobulk.drop_obs
   Pseudobulk.drop_var
   Pseudobulk.drop_cell_types
   Pseudobulk.rename_obs
   Pseudobulk.rename_var
   Pseudobulk.rename_cell_types
   Pseudobulk.cast_X
   Pseudobulk.cast_obs
   Pseudobulk.cast_var
   Pseudobulk.join_obs
   Pseudobulk.join_var
   Pseudobulk.subsample_obs
   Pseudobulk.subsample_var
   Pseudobulk.split_by_cell_type
   Pseudobulk.concat_obs
   Pseudobulk.concat_var

Transformation
--------------

.. autosummary::
   :nosignatures:

   Pseudobulk.copy
   Pseudobulk.to_df
   Pseudobulk.map_X
   Pseudobulk.map_obs
   Pseudobulk.map_var

Analysis
--------

.. autosummary::
   :nosignatures:

   Pseudobulk.qc
   Pseudobulk.library_size
   Pseudobulk.CPM
   Pseudobulk.log_CPM
   Pseudobulk.regress_out
   Pseudobulk.DE

Utility
-------

.. autosummary::
   :nosignatures:

   Pseudobulk.peek_obs
   Pseudobulk.peek_var
   Pseudobulk.pipe
   Pseudobulk.pipe_X
   Pseudobulk.pipe_obs
   Pseudobulk.pipe_var
