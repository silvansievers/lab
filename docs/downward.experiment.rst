:mod:`downward.experiment` --- Fast Downward experiment
=======================================================
.. autoclass:: downward.experiment.FastDownwardExperiment
   :members: add_algorithm, add_suite

Built-in parsers
----------------

The following constants are paths to built-in parsers that can be passed
to :meth:`exp.add_parser() <lab.experiment.Experiment.add_parser>`. See
:class:`FastDownwardExperiment
<downward.experiment.FastDownwardExperiment>` for examples.

Mandatory built-in parsers
~~~~~~~~~~~~~~~~~~~~~~~~~~

These two parsers have to be part of every Fast Downward experiment and
``LAB_DRIVER_PARSER`` has to be added before ``EXITCODE_PARSER``.

.. autoattribute:: lab.experiment.Experiment.LAB_DRIVER_PARSER
   :annotation:

.. autoattribute:: downward.experiment.FastDownwardExperiment.EXITCODE_PARSER
   :annotation:

Optional built-in parsers
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoattribute:: downward.experiment.FastDownwardExperiment.TRANSLATOR_PARSER
   :annotation:

.. autoattribute:: downward.experiment.FastDownwardExperiment.SINGLE_SEARCH_PARSER
   :annotation:

.. autoattribute:: downward.experiment.FastDownwardExperiment.ANYTIME_SEARCH_PARSER
   :annotation:
