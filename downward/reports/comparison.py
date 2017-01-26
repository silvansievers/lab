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
        for index, algo_pair in enumerate(algorithm_pairs):
            assert type(algo_pair) is tuple
            for algo in algo_pair[:2]:
                algos.add(algo)
        kwargs['filter_algorithm'] = algos
        self.algorithm_pairs = algorithm_pairs
        AbsoluteReport.__init__(self, **kwargs)
        self.algopair_domain_attribute_better = defaultdict(int)
        self.algopair_domain_attribute_worse = defaultdict(int)

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
                table.add_cell(domain, table._get_diff_col_name(algo_pair), diff)
                table.add_cell(domain, table._get_better_col_name(algo_pair),
                               self.algopair_domain_attribute_better[algo_pair, domain, attribute])
                table.add_cell(domain, table._get_worse_col_name(algo_pair),
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
                table.add_cell(problem, table._get_diff_col_name(algo_pair), diff)
                table.add_cell(problem, table._get_better_col_name(algo_pair), better)
                table.add_cell(problem, table._get_worse_col_name(algo_pair), worse)

        return table

    def _get_empty_table(self, attribute=None, title=None):
        """Return an empty table."""
        if title is None:
            assert attribute is not None
            title = attribute
            if self.output_format == 'tex':
                title = title.capitalize().replace('_', ' ')

        if attribute is not None and self.attribute_is_numeric(attribute):
            # Decide whether we want to highlight minima or maxima.
            kwargs = dict(
                min_wins=attribute.min_wins,
                colored=self.colored and attribute.min_wins is not None,
                digits=attribute.digits)
        else:
            # Do not highlight anything.
            kwargs = {}
        table = ComparisonTable(title, self.algorithm_pairs, **kwargs)
        link = '#%s' % title
        formatter = reports.CellFormatter(link=link)
        table.cell_formatters[table.header_row][table.header_column] = formatter
        return table

class ComparisonTable(reports.Table):
    def __init__(self, title, algorithm_pairs, min_wins=None, colored=False, digits=2):
        reports.Table.__init__(self, title, min_wins, colored, digits)
        self.algorithm_pairs = algorithm_pairs
        self.algopair_columnname = {}
        columns = []
        for index, algo_pair in enumerate(algorithm_pairs):
            assert type(algo_pair) is tuple
            self.algopair_columnname[algo_pair] = (
                'Diff%d' % index, 'Better%d' % index, 'Worse%d' % index)
            columns.extend(self.get_col_names_for_algorithm_pair(algo_pair))
        self.set_column_order(columns)

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

    def get_col_names_for_algorithm_pair(self, algo_pair):
        return [
            algo_pair[0],
            algo_pair[1],
            self._get_diff_col_name(algo_pair),
            self._get_better_col_name(algo_pair),
            self._get_worse_col_name(algo_pair)
        ]

    def _get_values_for_col_name(self, col_name):
        values = []
        for _row_name in self.row_names:
            value = self[_row_name].get(col_name)
            if value is not None:
                values.append(value)
        return values

    def get_summary_rows(self):
        """
        Returns a dictionary mapping names of summary rows to dictionaries
        mapping column names to values.
        """
        summary_rows = {}
        for row_name in self.summary_row_order:
            func = self.summary_funcs[row_name]
            summary_row = {}
            for algo_pair in self.algorithm_pairs:
                col_names = self.get_col_names_for_algorithm_pair(algo_pair)
                values_col1 = self._get_values_for_col_name(col_names[0])
                values_col2 = self._get_values_for_col_name(col_names[1])
                values_col4 = self._get_values_for_col_name(col_names[3])
                values_col5 = self._get_values_for_col_name(col_names[4])
                col1_aggregated_value = func(values_col1)
                col2_aggregated_value = func(values_col2)
                # Algo1 and Algo2 columns: use aggregation function
                summary_row[col_names[0]] = col1_aggregated_value
                summary_row[col_names[1]] = col2_aggregated_value
                # Diff column: difference = algo2 - algo1
                summary_row[col_names[2]] = col2_aggregated_value - col1_aggregated_value
                # Better and Worse columns: use sum over the column
                summary_row[col_names[3]] = sum(values_col4)
                summary_row[col_names[4]] = sum(values_col5)

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

        min_wins = self.get_min_wins(row_name)

        for algo_pair in self.algorithm_pairs:
            col_names = self.get_col_names_for_algorithm_pair(algo_pair)
            values = [row[col_names[i]] for i in range(5)]
            colors = ['gray' for i in range(5)]
            bolds = [False for i in range(5)]
            if min_wins is not None and values[0] is not None and values[1] is not None:
                # Highlight the value and the diff columns if the difference is non-zero
                if (min_wins and values[1] < values[0]) or (not min_wins and values[1] > values[0]):
                    if self.colored:
                        colors[0] = 'blue'
                        colors[1] = 'blue'
                        colors[2] = 'green'
                    else:
                        #bolds[1] = True
                        bolds[2] = True

                if (min_wins and values[0] < values[1]) or (not min_wins and values[0] > values[1]):
                    if self.colored:
                        colors[0] = 'blue'
                        colors[1] = 'blue'
                        colors[2] = 'red'
                    #else:
                        #bolds[0] = True

                # Always highlight better/worse columns (even if the difference
                # coincidentally is 0)
                if self.colored:
                    if values[3]: # better column: color green if non-zero
                        colors[3] = 'green'
                    if values[4]: # worse column: color red if non-zero
                        colors[4] = 'red'
                else:
                    if values[3]: # better column: color green if non-zero
                        bolds[3] = True
                    if values[4]: # worse column: color red if non-zero
                        bolds[4] = True

            for index, col_name in enumerate(col_names):
                row[col_name] = self._format_cell(
                    row_name, col_name, values[index],
                    color=colors[index], bold=bolds[index])
