# -*- coding: utf-8 -*-
#
# downward uses the lab package to conduct experiments with the
# Fast Downward planning system.
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

from collections import defaultdict
import logging

from lab import reports
from lab import tools
from lab.reports.markup import Document, ESCAPE_WORDBREAK

from downward.reports.absolute import AbsoluteReport

class ComparisonReport(AbsoluteReport):
    def __init__(self, algorithm_pairs, **kwargs):
        if 'filter_algorithm' in kwargs:
            logging.critical(
                'ComparisonReport doesn\'t support "filter_algorithm". '
                'Use "algorithm_pairs" to select and order algorithms.')
        if not algorithm_pairs:
            logging.critical(
                'ComparisonReport requires a list of algorithm pairs: '
                '"algorithm_pairs"')

        algos = set()
        self.algopair_columnname = {}
        for index, algo_pair in enumerate(algorithm_pairs):
            assert type(algo_pair) is tuple
            self.algopair_columnname[algo_pair] = (
                'Diff%d' % index, 'Better%d' % index, 'Worse%d' % index)
            for algo in algo_pair[:2]:
                algos.add(algo)
        kwargs['filter_algorithm'] = algos
        self.algorithm_pairs = algorithm_pairs
        AbsoluteReport.__init__(self, **kwargs)
        self.algopair_domain_attribute_better = defaultdict(int)
        self.algopair_domain_attribute_worse = defaultdict(int)

    # TODO: can we use the full name diff-%s-%s and still display a custom
    # value in the table header, e.g. only "diff"?
    def _get_diff_col_name(self, algo_pair):
        #return 'Diff-%s-%s' % (algo_pair[0], algo_pair[1])
        return self.algopair_columnname[algo_pair][0]

    def _get_better_col_name(self, algo_pair):
        #return '%s-better-than-%s' % (algo_pair[1], algo_pair[0])
        return self.algopair_columnname[algo_pair][1]

    def _get_worse_col_name(self, algo_pair):
        #return '%s-worse-than-%s' % (algo_pair[1], algo_pair[0])
        return self.algopair_columnname[algo_pair][2]

    def get_markup(self):
        sections = []
        toc_lines = []

        warnings = self._get_warnings_table()
        if warnings:
            toc_lines.append('- **[''Unexplained Errors'' #unexplained-errors]**')
            sections.append(('unexplained-errors', warnings))

        toc_lines.append('- **[Info #info]**')
        sections.append(('info', self._get_general_info()))

        # Index of summary section.
        summary_index = len(sections)

        # Build a table containing summary functions of all other tables.
        # The actual section is added at position summary_index after creating
        # all other tables.
        summary = self._get_empty_table(title='Summary')
        summary.colored = self.colored
        toc_lines.append('- **[Summary #summary]**')

        for attribute in self.attributes:
            logging.info('Creating table(s) for %s' % attribute)
            tables = []
            ### Silvan: modified compared to absolute report
            # First create the domain tables...
            suite_table_index = len(tables)
            for domain in sorted(self.domains.keys()):
                tables.append((domain, self._get_table(attribute, domain)))

            # ... and then the suite table, because it needs values computed
            # during the creation of the domain tables. Insert the table before
            # the domain tables.
            if self.attribute_is_numeric(attribute):
                domain_table = self._get_table(attribute)
                tables.insert(suite_table_index, ('', domain_table))
                reports.extract_summary_rows(
                    domain_table, summary, link='#' + attribute)
            else:
                tables.insert(suite_table_index, (
                    '',
                    'Domain-wise reports only support numeric '
                    'attributes, but %s has type %s.' %
                    (attribute, self._all_attributes[attribute].__name__)))
            ### end modifications

            parts = []
            toc_line = []
            for (domain, table) in tables:
                if domain:
                    assert table
                    toc_line.append("[''%(domain)s'' #%(attribute)s-%(domain)s]" %
                                    locals())
                    parts.append('== %(domain)s ==[%(attribute)s-%(domain)s]\n'
                                 '%(table)s\n' % locals())
                else:
                    if table:
                        parts.append('%(table)s\n' % locals())
                    else:
                        parts.append('No task was found where all algorithms '
                                     'have a value for "%s". Therefore no '
                                     'domain-wise table can be generated.\n' %
                                     attribute)

            toc_lines.append("- **[''%s'' #%s]**" % (attribute, attribute))
            toc_lines.append('  - ' + ' '.join(toc_line))
            sections.append((attribute, '\n'.join(parts)))

        # Add summary before main content. This is done after creating the main content
        # because the summary table is extracted from all other tables.
        sections.insert(summary_index, ('summary', summary))

        toc = '\n'.join(toc_lines)

        content = '\n'.join('= %s =[%s]\n\n%s' % (attr, attr, section)
                            for (attr, section) in sections)
        return '%s\n\n\n%s' % (toc, content)

    def _get_suite_table(self, attribute):
        assert self.attribute_is_numeric(attribute), attribute
        table = self._get_empty_table(attribute)
        self._add_summary_functions(table, attribute)
        # The first group function is used for aggregation.
        func_name, func = self._get_group_functions(attribute)[0]
        num_probs = 0
        self._add_table_info(attribute, func_name, table)
        domain_algo_values = defaultdict(list)
        for (domain, problem), runs in self.problem_runs.items():
                if (not attribute.absolute and
                        any(run.get(attribute) is None for run in runs)):
                    continue
                num_probs += 1
                for run in runs:
                    value = run.get(attribute)
                    if value is not None:
                        domain_algo_values[(domain, run['algorithm'])].append(value)

        # If the attribute is absolute (e.g. coverage) we may have
        # added problems for which not all algorithms have a value. Therefore, we
        # can only print the number of instances (in brackets after the domain
        # name) if that number is the same for all algorithms. If not all algorithms
        # have values for the same number of problems, we write the full list of
        # different problem numbers.
        num_values_lists = defaultdict(list)
        for domain in self.domains:
            for algo in self.algorithms:
                values = domain_algo_values.get((domain, algo), [])
                num_values_lists[domain].append(str(len(values)))
        for domain, num_values_list in num_values_lists.items():
            if len(set(num_values_list)) == 1:
                count = num_values_list[0]
            else:
                count = ','.join(num_values_list)
            link = None
            if self.use_domain_links:
                link = '#%s-%s' % (attribute, domain)
            formatter = reports.CellFormatter(link=link, count=count)
            table.cell_formatters[domain][table.header_column] = formatter

        ### Silvan: modified compared to absolute report
        domain_algo_aggregatedresult = {}
        for (domain, algo), values in domain_algo_values.items():
            aggreated_result = func(values)
            table.add_cell(domain, algo, aggreated_result)
            domain_algo_aggregatedresult[domain, algo] = aggreated_result

        table.num_values = num_probs

        for algo_pair in self.algorithm_pairs:
            for domain in self.domains:
                algo2_value = domain_algo_aggregatedresult.get((domain, algo_pair[1]), None)
                algo1_value = domain_algo_aggregatedresult.get((domain, algo_pair[0]), None)
                if algo2_value is None or algo1_value is None:
                    continue
                diff = algo2_value - algo1_value
                table.add_cell(domain, self._get_diff_col_name(algo_pair), diff)
                table.add_cell(domain, self._get_better_col_name(algo_pair),
                               self.algopair_domain_attribute_better[algo_pair, domain, attribute])
                table.add_cell(domain, self._get_worse_col_name(algo_pair),
                               self.algopair_domain_attribute_worse[algo_pair, domain, attribute])
        ### end modifications

        return table

    def _get_domain_table(self, attribute, domain):
        table = self._get_empty_table(attribute)

        min_wins = attribute.min_wins
        for algo_pair in self.algorithm_pairs:
            algo1 = algo_pair[0]
            algo2 = algo_pair[1]
            for problem in self.domains[domain]:
                algo1_value = self.runs[domain, problem, algo1].get(attribute)
                algo2_value = self.runs[domain, problem, algo2].get(attribute)
                try:
                    diff = float(algo2_value) - float(algo1_value)
                except (ValueError, TypeError, KeyError):
                    diff = None
                better = 0
                worse = 0
                if diff is not None:
                    if (diff < 0 and min_wins) or (diff > 0 and not min_wins): # algo2 is better
                        better = 1
                        self.algopair_domain_attribute_better[algo_pair, domain, attribute] += 1
                    elif (diff > 0 and min_wins) or (diff < 0 and not min_wins): # algo2 is worse
                        worse = 1
                        self.algopair_domain_attribute_worse[algo_pair, domain, attribute] += 1

                table.add_cell(problem, algo1, algo1_value)
                table.add_cell(problem, algo2, algo2_value)
                table.add_cell(problem, self._get_diff_col_name(algo_pair), diff)
                table.add_cell(problem, self._get_better_col_name(algo_pair), better)
                table.add_cell(problem, self._get_worse_col_name(algo_pair), worse)

        return table

    def _get_empty_table(self, attribute=None, title=None):
        """Return an empty table."""
        if title is None:
            assert attribute is not None
            title = attribute
            if self.output_format == 'tex':
                title = title.capitalize().replace('_', ' ')

        columns = []
        for algo_pair in self.algorithm_pairs:
            columns.append(algo_pair[0])
            columns.append(algo_pair[1])
            columns.append(self._get_diff_col_name(algo_pair))
            columns.append(self._get_better_col_name(algo_pair))
            columns.append(self._get_worse_col_name(algo_pair))

        if attribute is not None and self.attribute_is_numeric(attribute):
            # Decide whether we want to highlight minima or maxima.
            kwargs = dict(
                min_wins=attribute.min_wins,
                colored=self.colored and attribute.min_wins is not None,
                digits=attribute.digits)
        else:
            # Do not highlight anything.
            kwargs = {}
        table = ComparisonTable(title=title, **kwargs)
        table.set_column_order(columns)
        for algo_pair in self.algorithm_pairs:
            table.set_column_color(self._get_diff_col_name(algo_pair), 'diff')
            table.set_column_color(self._get_better_col_name(algo_pair), 'better')
            table.set_column_color(self._get_worse_col_name(algo_pair), 'worse')
        link = '#%s' % title
        formatter = reports.CellFormatter(link=link)
        table.cell_formatters[table.header_row][table.header_column] = formatter
        return table

class ComparisonTable(reports.Table):
    def __init__(self, title='', min_wins=None, colored=False, digits=2):
        reports.Table.__init__(self, title, min_wins, colored, digits)
        self.column_color_type = {}

    def set_column_color(self, col_name, color):
        self.column_color_type[col_name] = color

    def get_summary_rows(self):
        """
        Returns a dictionary mapping names of summary rows to dictionaries
        mapping column names to values.
        """
        summary_rows = {}
        for row_name in self.summary_row_order:
            func = self.summary_funcs[row_name]
            summary_row = {}
            second_to_previous_col_aggregated_value = -1
            previous_col_aggregated_value = -1
            col_counter = 0
            for col_name in self.col_names:
                values = []
                for _row_name in self.row_names:
                    value = self[_row_name].get(col_name)
                    if value is not None:
                        values.append(value)
                if self.column_color_type.get(col_name, '') == 'diff':
                    # For diff columns, compute the diff of the previously stored
                    # aggregated values of the first and second column.
                    assert second_to_previous_col_aggregated_value != -1
                    assert previous_col_aggregated_value != -1
                    summary_row[col_name] = previous_col_aggregated_value - second_to_previous_col_aggregated_value
                else:
                    if values:
                        # Always use sum to aggregate better/worse tables
                        if self.column_color_type.get(col_name, '') in ['better', 'worse']:
                            summary_row[col_name] = sum(values)
                        else:
                            summary_row[col_name] = func(values)
                    else:
                        summary_row[col_name] = None
                # TODO: this assumes a fixed number of 5 columns for every comparison
                if col_counter % 5 == 0:
                    second_to_previous_col_aggregated_value = summary_row[col_name]
                if col_counter % 5 == 1:
                    previous_col_aggregated_value = summary_row[col_name]
                col_counter += 1

            summary_row[self.header_column] = row_name
            summary_rows[row_name] = summary_row
            formatter = reports.CellFormatter(bold=True, count=self.num_values)
            self.cell_formatters[row_name][self.header_column] = formatter
        return summary_rows

    def _format_row(self, row_name, row):
        """Format all entries in **row** (in place)."""
        if row_name == self.header_row:
            for col_name, value in row.items():
                # Allow breaking after underlines.
                value = value.replace('_', '_' + ESCAPE_WORDBREAK)
                # Right-align headers (except the left-most one).
                if col_name != self.header_column:
                    value = ' ' + value
                row[col_name] = value
            return

        # Get the slice of the row that should be formated (i.e. the data columns).
        # Note that there might be other columns (e.g. added by dynamic data
        # modules) that should not be formated.
        row_slice = dict((col_name, row.get(col_name))
                         for col_name in self.col_names if self.column_color_type.get(col_name, None) is None)

        min_value, max_value = tools.get_min_max(row_slice.values())

        min_wins = self.get_min_wins(row_name)
        highlight = min_wins is not None
        colors = tools.get_colors(row_slice, min_wins) if self.colored else None

        def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
            return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

        for col_name, value in row.items():
            color = None
            bold = False
            # Format data columns
            if col_name in row_slice:
                if self.colored:
                    color = tools.rgb_fractions_to_html_color(*colors[col_name])
                elif highlight and value is not None and (
                        (is_close(value, min_value) and min_wins) or
                        (is_close(value, max_value) and not min_wins)):
                    bold = True
            # Decide on the color of diff/better/worse columns based on value
            elif self.column_color_type.get(col_name, None):
                color_type = self.column_color_type[col_name]
                if color_type == 'diff':
                    if value is None or is_close(value, 0):
                        color = 'gray'
                    elif (value < 0 and min_wins) or (value > 0 and not min_wins):
                        color = 'green'
                    elif (value > 0 and min_wins) or (value < 0 and not min_wins):
                        color = 'red'
                    else:
                        assert False
                elif color_type == 'better':
                    if value == 0:
                        color = 'gray'
                    else:
                        color = 'green'
                elif color_type == 'worse':
                    if value == 0:
                        color = 'gray'
                    else:
                        color = 'red'
                else:
                    assert False
            row[col_name] = self._format_cell(row_name, col_name, value,
                                              color=color, bold=bold)


