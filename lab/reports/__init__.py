# -*- coding: utf-8 -*-
#
# lab is a Python API for running and evaluating algorithms.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Module that permits generating reports by reading properties files
"""

from __future__ import with_statement, division

import collections
from collections import defaultdict
import fnmatch
import logging
import math
import numbers
import os

from lab import tools
from lab.external import txt2tags
from lab.reports.markup import Document, ESCAPE_WORDBREAK


@tools.remove_none_values
def prod(values):
    """Computes the product of a list of numbers.

    >>> print prod([2, 3, 7])
    42.0
    """
    prod = 1
    for value in values:
        prod *= value
    return prod


@tools.remove_none_values
def avg(values):
    """Compute the arithmetic mean of a list of numbers.

    >>> avg([20, 30, 70])
    40.0
    """
    return round(math.fsum(values) / len(values), 4)


@tools.remove_none_values
def gm(values):
    """Compute the geometric mean of a list of numbers.

    >>> gm([2, 8])
    4.0
    """
    values = [val if val > 0 else 0.1 for val in values]
    exp = 1.0 / len(values)
    return round(prod([val ** exp for val in values]), 4)


@tools.remove_none_values
def minimum(values):
    """Filter out None values and return the minimum.

    If there are only None values, return None.
    """
    return min(values)


@tools.remove_none_values
def maximum(values):
    """Filter out None values and return the maximum.

    If there are only None values, return None.
    """
    return max(values)


@tools.remove_none_values
def stddev(values):
    """Compute the standard deviation of a list of numbers.

    >>> stddev([2, 4, 4, 4, 5, 5, 7, 9])
    2.0
    """
    n = len(values)
    mu = avg(values)
    return math.sqrt((sum((v - mu) ** 2 for v in values) / n))


def function_name(f):
    names = {'avg': 'average', 'gm': 'geometric mean'}
    return names.get(f.__name__, f.__name__)


class Attribute(str):
    """A string subclass for attributes in reports."""
    def __new__(cls, name, **kwargs):
        return str.__new__(cls, name)

    def __init__(self, name, absolute=False, min_wins=True, functions=sum):
        """Use this class if your **own** attribute needs a non-default value for:

        * *absolute*: If False, only include tasks for which all task runs have
          values in a domain-wise table (e.g. ``coverage`` is absolute, whereas
          ``expansions`` is not, because we can't compare configurations A and B
          for task X if B has no value for ``expansions``).
        * *min_wins*: Set to True if a smaller value for this attribute is
          better and to False otherwise (e.g. for ``coverage`` *min_wins* is
          False, whereas it is True for ``expansions``).
        * *functions*: Set the function or functions that are used to group values
          of multiple runs for this attribute. The first entry is used to aggregate
          values for domain-wise reports (e.g. for ``coverage`` this is
          :py:func:`sum`, whereas ``expansions`` uses :py:func:`gm`). This can be a
          single function or a list of functions and defaults to :py:func:`sum`.

        The ``downward`` package automatically uses appropriate settings for
        most attributes. ::

            from lab.reports import minimum, maximum
            avg_h = Attribute('average_h', min_wins=False,
                              functions=[sum, minimum, maximum])
            abstraction_done = Attribute('abstraction_done', absolute=True)

            Report(attributes=[avg_g, abstraction_done, 'coverage', 'expansions'])

        """
        self.absolute = absolute
        self.min_wins = min_wins
        if not isinstance(functions, collections.Iterable):
            functions = [functions]
        self.functions = functions

    def copy(self, name):
        return Attribute(name, absolute=self.absolute, min_wins=self.min_wins,
                         functions=self.functions)


class Report(object):
    """
    Base class for all reports.
    """
    def __init__(self, attributes=None, format='html', filter=None, **kwargs):
        """
        *attributes* is a list of the attributes you want to include in your
        report. If omitted, use all found numerical attributes. Globbing
        characters * and ? are allowed. Example: ::

            Report(attributes=['translator_time_*'])

        When a report is made, both the available and the selected attributes
        are printed on the commandline.

        *format* can be one of e.g. html, tex, wiki (MediaWiki),
        gwiki (Google Code Wiki), doku (DokuWiki), pmw (PmWiki),
        moin (MoinMoin), txt (Plain text) and art (ASCII art). Subclasses may
        allow additional formats.

        If given, *filter* must be a function or a list of functions that
        are passed a dictionary of a run's keys and values and return
        True or False. Depending on the returned value, the run is included
        or excluded from the report.
        Alternatively, the function can return a dictionary that will overwrite
        the old run's dictionary for the report.

        Filters for properties can be given in shorter form without defining a function
        To include only runs where property p has value v, use *filter_p=v*.
        To include only runs where property p has value v1, v2 or v3, use
        *filter_p=[v1, v2, v3]*.

        Examples:

        Include only *coverage* and *expansions* in the report and write a
        LaTeX file at ``myreport.tex``::

            report = Report(attributes=['coverage', 'expansions'], format='tex')
            report(path_to_eval_dir, 'myreport.tex')

        Only include successful runs in the report::

            report = Report(filter_coverage=1)
            report(path_to_eval_dir, 'myreport.html')

        Only include runs in the report where the time score is better than the
        memory score::

            def better_time_than_memory_score(run):
                return run['score_search_time'] > run['score_memory']
            report = Report(filter=better_time_than_memory_score)
            report(path_to_eval_dir, 'myreport.html')

        Filter function that filters and renames configs with additional sorting::

            def rename_configs(run):
                config = run['config'].replace('WORK-', '')
                paper_names = {'lama11': 'LAMA 2011', 'fdss_sat1': 'FDSS 1',
                               'fdss_sat2': 'FDSS 2'}
                run['config'] = paper_names.get(config, 'unknown')
                return run

            # We want LAMA 2011 to be the leftmost column.
            # Filters defined with key word arguments are evaluated last,
            # so we use the updated config names here.
            configs = ['LAMA 2011', 'FDSS 1', 'FDSS 2']
            Report(filter=rename_configs, filter_config=configs)

        Filter function that only allows runs with a timeout in one of two domains::

            report = Report(attributes=['coverage'],
                            filter_domain=['blocks', 'barman'],
                            filter_search_timeout=1)
        """
        if isinstance(attributes, basestring):
            attributes = [attributes]
        self.attributes = attributes or []
        assert format in txt2tags.TARGETS + ['eps', 'pdf', 'pgf', 'png', 'py']
        self.output_format = format
        self.toc = True
        self.run_filter = tools.RunFilter(filter, **kwargs)

    def __call__(self, eval_dir, outfile):
        """Make the report.

        *eval_dir* must be a path to an evaluation directory containing a
        ``properties`` file.

        The report will be written to *outfile*.
        """
        if not eval_dir.endswith('-eval'):
            logging.info('The source directory does not end with "-eval". '
                         'Are you sure this is an evaluation directory?')
        self.eval_dir = os.path.abspath(eval_dir)
        # It would be nice if we could infer "format" from "outfile", but the
        # former is needed before the latter is available.
        # Also we cannot add the extension ".format" to "outfile" in case it's
        # missing, because "outfile" might be a directory.
        self.outfile = os.path.abspath(outfile)

        # Map from attribute to type.
        self._all_attributes = {}
        self._load_data()
        self._apply_filter()
        self._scan_data()

        # Turn string attributes into instances of Attribute.
        self.attributes = [self._prepare_attribute(attr)
                           for attr in self.attributes]

        # Expand glob characters.
        self.attributes = self._glob_attributes(self.attributes)

        if not self.attributes:
            logging.info('Available attributes: %s' % ', '.join(self.all_attributes))
            logging.info('Using all numerical attributes.')
            self.attributes = self._get_numerical_attributes()

        self.attributes = sorted(set(self.attributes))
        self.write()

    def _prepare_attribute(self, attr):
        if isinstance(attr, Attribute):
            return attr
        return Attribute(attr)

    def _glob_attributes(self, attributes):
        expanded_attrs = []
        for attr in attributes:
            # Attribute without wildcards. Filtering would reset its options.
            if attr in self.all_attributes:
                expanded_attrs.append(attr)
                continue
            matches = fnmatch.filter(self.all_attributes, attr)
            if not matches:
                logging.warning('Attribute %s is not present in the dataset.' %
                                attr)
            # Use the attribute options from the pattern for all matches, but
            # don't try to guess options for attributes that appear in the list.
            expanded_attrs.extend([attr.copy(match) for match in matches
                                   if not match in attributes])
        if attributes and not expanded_attrs:
            logging.critical('No attributes match your patterns.')
        return expanded_attrs

    @property
    def all_attributes(self):
        return sorted(self._all_attributes.keys())

    def _get_numerical_attributes(self):
        return [attr for attr in self._all_attributes.keys()
                if self.attribute_is_numeric(attr)]

    def attribute_is_numeric(self, attribute):
        """Return true if the values for *attribute* are ints or floats.

        If the attribute is None in all runs it may be numeric.

        """
        return (self._all_attributes[attribute] is None or
                issubclass(self._all_attributes[attribute], numbers.Number))

    def get_markup(self):
        """
        If ``get_text()`` is not overwritten, this method can be overwritten in
        subclasses that want to return the report as
        `txt2tags <http://txt2tags.org/>`_ markup. The default ``get_text()``
        method converts that markup into *format*.
        """
        table = Table()
        for run_id, run in self.props.items():
            row = {}
            for key, value in run.items():
                if not key in self.attributes:
                    continue
                if isinstance(value, (list, tuple)):
                    key = '-'.join([str(item) for item in value])
                row[key] = value
            table.add_row(run_id, row)
        return str(table)

    def get_text(self):
        """
        This method should be overwritten in subclasses that want to produce
        e.g. HTML or LaTeX markup or programming code directly instead of
        creating `txt2tags <http://txt2tags.org/>`_ markup.
        """
        name, ext = os.path.splitext(os.path.basename(self.outfile))
        doc = Document(title=name)
        markup = self.get_markup()

        if not markup:
            markup = ('No tables were generated. '
                      'This happens when no significant changes occured or '
                      'if for all attributes and all problems never all '
                      'configs had a value for this attribute in a '
                      'domain-wise report.')

        doc.add_text(markup)
        if len(markup) < 100000:
            print 'REPORT MARKUP:\n'
            print doc
        return doc.render(self.output_format, {'toc': self.toc})

    def write(self):
        """
        Overwrite this method if you want to write the report directly. You
        should write the report to *self.outfile*.
        """
        content = self.get_text()
        tools.makedirs(os.path.dirname(self.outfile))
        with open(self.outfile, 'w') as file:
            file.write(content)
            logging.info('Wrote file://%s' % self.outfile)

    def _get_type(self, attribute):
        for run_id, run in self.props.iteritems():
            val = run.get(attribute)
            if val is not None:
                return type(val)
        # Attribute is None in all runs.
        return None

    def _get_type_map(self, attributes):
        return dict((self._prepare_attribute(attr), self._get_type(attr))
                    for attr in attributes)

    def _scan_data(self):
        attributes = set()
        for run_id, run in self.props.items():
            attributes |= set(run.keys())
        self._all_attributes = self._get_type_map(attributes)

    def _load_data(self):
        props_file = os.path.join(self.eval_dir, 'properties')
        if not os.path.exists(props_file):
            logging.critical('Properties file not found at %s' % props_file)

        logging.info('Reading properties file')
        self.props = tools.Properties(filename=props_file)
        logging.info('Reading properties file finished')
        if not self.props:
            logging.critical('properties file in evaluation dir is empty.')

    def _apply_filter(self):
        self.props = self.run_filter.apply(self.props)
        if not self.props:
            logging.critical('All runs have been filtered -> Nothing to report.')


class CellFormatter(object):
    """Formating information for one cell in a table."""
    def __init__(self, bold=False, count=None, link=None):
        self.bold = bold
        self.count = count
        self.link = link

    def format_value(self, value):
        result = str(value)
        if self.link:
            result = '[""%s"" %s]' % (result, self.link)
        if self.count:
            result = '%s (%s)' % (result, self.count)
        if self.bold:
            result = '**%s**' % result
        return result


class Table(collections.defaultdict):
    def __init__(self, title='', min_wins=None, colored=False):
        """
        The *Table* class can be useful for `Report` subclasses that want to
        return a table as txt2tags markup. It is realized as a dictionary of
        dictionaries mapping row names to colum names to cell values. To obtain
        the markup from a table, use the ``str()`` function.

        *title* will be printed in the top left cell.

        *min_wins* can be either None, True or False. If it is True (False),
        the cell with the lowest (highest) value in each row will be
        highlighted.

        If *colored* is True, the values of each row will be given colors from a
        colormap.

        >>> t = Table(title='expansions')
        >>> t.add_cell('prob1', 'cfg1', 10)
        >>> t.add_cell('prob1', 'cfg2', 20)
        >>> t.add_row('prob2', {'cfg1': 15, 'cfg2': 25})
        >>> print t
        || expansions |  cfg1 |  cfg2 |
         | prob1      |    10 |    20 |
         | prob2      |    15 |    25 |
        >>> t.row_names
        ['prob1', 'prob2']
        >>> t.col_names
        ['cfg1', 'cfg2']
        >>> t.get_row('prob2')
        [15, 25]
        >>> t.get_columns() == {'cfg1': [10, 15], 'cfg2': [20, 25]}
        True
        >>> t.add_summary_function('SUM', sum)
        >>> print t
        || expansions |  cfg1 |  cfg2 |
         | prob1      |    10 |    20 |
         | prob2      |    15 |    25 |
         | **SUM**    |    25 |    45 |
        >>> t.set_column_order(['cfg2', 'cfg1'])
        >>> print t
        || expansions |  cfg2 |  cfg1 |
         | prob1      |    20 |    10 |
         | prob2      |    25 |    15 |
         | **SUM**    |    45 |    25 |
        """
        collections.defaultdict.__init__(self, dict)

        self.title = title
        self.min_wins = min_wins
        self.row_min_wins = {}
        self.colored = colored

        self.summary_funcs = {}
        self.info = []
        self.num_values = None
        self.dynamic_data_modules = []

        self._cols = None

        # For printing.
        self.header_row = 'column names (never printed)'
        self.header_column = 'row names (never printed)'
        self.cell_formatters = collections.defaultdict(dict)
        self.col_size = None
        self.column_order = None
        self.summary_row_order = []

    def add_cell(self, row, col, value):
        """Set Table[row][col] = value."""
        self[row][col] = value
        self._cols = None

    def add_row(self, row_name, row):
        """Add a new data row called *row_name* to the table.

        *row* must be a mapping from column names to values.
        """
        self[row_name] = row
        self._cols = None

    def add_col(self, col_name, col):
        """Add a new data column called *col_name* to the table.

        *col* must be a mapping from row names to values.
        """
        for row_name, value in col.items():
            self[row_name][col_name] = value
        self._cols = None

    @property
    def row_names(self):
        """Return all data row names in sorted order."""
        return tools.natural_sort(self.keys())

    @property
    def col_names(self):
        """Return all data column names in sorted order."""
        if self._cols:
            return self._cols
        col_names = set()
        for row in self.values():
            col_names |= set(row.keys())
        self._cols = []
        if self.column_order:
            # First use all elements for which we know an order.
            # All remaining elements will be sorted alphabetically.
            self._cols = [c for c in self.column_order if c in col_names]
            col_names -= set(self._cols)
        self._cols += tools.natural_sort(col_names)
        return self._cols

    def get_row(self, row_name):
        """Return a list of the values in *row*."""
        return [self[row_name].get(col_name, None) for col_name in self.col_names]

    def get_columns(self):
        """
        Return a mapping from column names to the list of values in that column.
        """
        values = defaultdict(list)
        for row_name in self.row_names:
            for col_name in self.col_names:
                values[col_name].append(self[row_name].get(col_name))
        return values

    def add_summary_function(self, name, func):
        """
        Add a bottom row with the values ``func(column_values)`` for each column.
        *func* can be e.g. ``sum``, ``reports.avg`` or ``reports.gm``.
        """
        self.summary_funcs[name] = func
        self.summary_row_order.append(name)

    def set_column_order(self, order):
        self.column_order = order
        self._cols = None

    def get_min_wins(self, row_name=None):
        """
        The table class can store information on whether higher or lower values are better
        for each row or globally. If no row specific setting for *row_name* is found, the
        global setting is returned.
        """
        return self.row_min_wins.get(row_name, self.min_wins)

    def get_summary_rows(self):
        """
        Returns a dictionary mapping names of summary rows to dictionaries
        mapping column names to values.
        """
        summary_rows = {}
        for row_name in self.summary_row_order:
            func = self.summary_funcs[row_name]
            summary_row = {}
            for col_name, column in self.get_columns().items():
                values = [val for val in column if val is not None]
                if values:
                    summary_row[col_name] = func(values)
                else:
                    summary_row[col_name] = None
            summary_row[self.header_column] = row_name
            summary_rows[row_name] = summary_row
            formatter = CellFormatter(bold=True, count=self.num_values)
            self.cell_formatters[row_name][self.header_column] = formatter
        return summary_rows

    def _get_printable_row_order(self):
        """
        Return a list of all rows (including non-data rows) in the order
        they should be printed.
        """
        row_order = [self.header_row]
        for row_name in self.row_names + self.summary_row_order:
            row_order.append(row_name)
        for module in self.dynamic_data_modules:
            row_order = module.modify_printable_row_order(self, row_order) or row_order
        return row_order

    def _get_printable_column_order(self):
        """
        Return a list of all columns (including non-data columns) in the order
        they should be printed.
        """
        col_order = [self.header_column]
        for col_name in self.col_names:
            col_order.append(col_name)
        for module in self.dynamic_data_modules:
            col_order = module.modify_printable_column_order(self, col_order) or col_order
        return col_order

    def _collect_cells(self):
        """
        Collect all cells that should be printed including table headers,
        row names, summary rows, etc. Returns a dictionary mapping row names
        to dictionaries mapping column names to values.
        """
        cells = collections.defaultdict(dict)
        cells[self.header_row][self.header_column] = self.title
        for col_name in self.col_names:
            cells[self.header_row][col_name] = str(col_name)
        # Add data rows and summary rows.
        for row_name, row in self.items() + self.get_summary_rows().items():
            cells[row_name][self.header_column] = str(row_name)
            for col_name in self.col_names:
                cells[row_name][col_name] = row.get(col_name)
        for dynamic_data_module in self.dynamic_data_modules:
            cells = dynamic_data_module.collect(self, cells) or cells
        return cells

    def _format(self, cells):
        """Format all entries in **cells** (in place)."""
        for row_name, row in cells.items():
            self._format_row(row_name, row)
        for dynamic_data_module in self.dynamic_data_modules:
            dynamic_data_module.format(self, cells)

    def _format_row(self, row_name, row):
        """Format all entries in **row** (in place)."""
        if row_name == self.header_row:
            for col_name, value in row.items():
                row[col_name] = value.replace('_', '_' + ESCAPE_WORDBREAK)
            return

        # Get the slice of the row that should be formated (i.e. the data columns).
        # Note that there might be other columns (e.g. added by dynamic data
        # modules) that should not be formated.
        row_slice = dict((col_name, row.get(col_name))
                         for col_name in self.col_names)

        min_value, max_value = tools.get_min_max(row_slice.values())

        min_wins = self.get_min_wins(row_name)
        highlight = min_wins is not None
        colors = tools.get_colors(row_slice, min_wins) if self.colored else None

        for col_name, value in row.items():
            color = None
            bold = False
            # Format data columns
            if col_name in row_slice:
                rounded_value = round(value, 2) if isinstance(value, float) else value
                if self.colored:
                    color = tools.rgb_fractions_to_html_color(*colors[col_name])
                elif highlight and (rounded_value == min_value and min_wins or
                                    rounded_value == max_value and not min_wins):
                    bold = True
            row[col_name] = self._format_cell(row_name, col_name, value,
                                              color=color, bold=bold)

    def _format_cell(self, row_name, col_name, value, color=None, bold=False):
        """
        Return the formatted value for a single cell in the table.
        *row_name* and *col_name* specify the position of the cell and *value* is the
        unformatted value.
        Floats are rounded to two decimal places and lists are quoted. The *color* to
        render the result in can be given as a string and setting *bold* to true
        renders the value in bold.

        If a custom formatter is specified for this cell, it is used instead of this
        default format.
        """
        formatter = self.cell_formatters.get(row_name, {}).get(col_name)
        if formatter:
            return formatter.format_value(value)

        justify_right = isinstance(value, (float, int))

        def format_value(value):
            if isinstance(value, float):
                return '%.2f' % value
            elif isinstance(value, (list, tuple)):
                # Escape brackets to avoid involuntarily creating txt2tags links.
                return "''[''%s'']''" % (', '.join(format_value(v) for v in value) or ' ')
            else:
                return str(value)

        value_text = format_value(value)

        if color is not None:
            value_text = '{%s|color:%s}' % (value_text, color)
        if bold:
            value_text = '**%s**' % value_text
        if justify_right:
            value_text = ' ' + value_text
        return value_text

    def _get_markup(self, cells):
        """
        Return a string cotaining all printable cells (see
        **_get_printable_column_order** and **_get_printable_row_order**)
        as correctly formatted markup.
        """
        parts = []
        for row_name in self._get_printable_row_order():
            if row_name == self.header_row:
                parts.append(self._get_header_markup(row_name, cells[row_name]))
            else:
                parts.append(self._get_row_markup(row_name, cells[row_name]))
        if self.info:
            parts.append(' '.join(self.info))
        return '\n'.join(parts)

    def _get_header_markup(self, row_name, row):
        """Return the txt2tags table markup for the headers."""
        return self._get_row_markup(row_name, row, template='|| %s |')

    def _get_row_markup(self, row_name, row, template=' | %s |'):
        """Return the txt2tags table markup for one row."""
        formatted_cells = []
        for col_name in self._get_printable_column_order():
            cell = self._get_cell_markup(row_name, col_name, row.get(col_name, ''))
            formatted_cells.append(cell)
        return template % ' | '.join(formatted_cells)

    def _get_cell_markup(self, row_name, col_name, value):
        """Values are already formatted."""
        # TODO: Remove this method?
        return str(value)

    def __str__(self):
        """Return the txt2tags markup for this table."""
        cells = self._collect_cells()
        self._format(cells)
        return self._get_markup(cells)


def extract_summary_rows(from_table, to_table, link=None):
    """
    Extract all summary rows of **from_table** and add them as data rows
    to **to_table**.
    """
    for name, row in from_table.get_summary_rows().items():
        row_name = '%s - %s' % (from_table.title, name)
        if link is not None:
            formatter = CellFormatter(link=link)
            to_table.cell_formatters[row_name][to_table.header_column] = formatter
        to_table.row_min_wins[row_name] = from_table.min_wins
        for col_name, value in row.items():
            if col_name == from_table.header_column:
                continue
            to_table.add_cell(row_name, col_name, value)


class DynamicDataModule(object):
    """Interface for modules that dynamically add or modify data in a table."""
    def collect(self, table, cells):
        """
        Called after the data collection in the table. Subclasses can
        add new values to **cells** or modify existing values.
        """
        return cells

    def format(self, table, formated_cells):
        """
        Called after the formatting in the table. Subclasses can
        (re-)format all values in **formated_cells**. Specifically all new
        values added by the **collect** method should be formatted.
        """
        pass

    def modify_printable_row_order(self, table, row_order):
        """
        Called after retrieving a row order in the table. Subclasses can
        modify the order or add new rows. Specifically all rows that were
        added by the **collect** method should be appended or
        inserted.
        """
        return row_order

    def modify_printable_column_order(self, table, column_order):
        """
        Called after retrieving a column order in the table. Subclassed can
        modify the order or add new columns. Specifically all columns that were
        values added by the **collect** method should be appended or
        inserted.
        """
        return column_order
