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

from downward.reports import PlanningReport


class AbsoluteReport(PlanningReport):
    """
    Write an absolute report about the attribute attribute, e.g.

        +------------+--------+--------+
        | expansions | hFF    | hCEA   |
        +============+========+========+
        | gripper    | 118    | 72     |
        +------------+--------+--------+
        | zenotravel | 21     | 17     |
        +------------+--------+--------+

    This report should be part of all your Fast Downward experiments as it
    automatically generates a table of unexplained errors, e.g. invalid solutions,
    unexpected timeouts and memory overflows. You should make sure that you
    check where the errors come from.
    """
    def __init__(self, resolution='combined', colored=True, **kwargs):
        """
        *resolution* must be one of "domain" or "problem" or "combined" (default).

        If *colored* is True, the values of each row will be given colors from a
        colormap. Only HTML reports can be colored currently.
        """
        PlanningReport.__init__(self, **kwargs)
        assert resolution in ['domain', 'problem', 'combined']
        self.resolution = resolution
        if colored and 'html' not in self.output_format:
            logging.info('Only HTML reports can be colored. Setting colored=False.')
            colored = False
        self.colored = colored
        self.toc = False

    def get_markup(self):
        sections = []
        toc_lines = []

        warnings = self._get_warnings_table()
        if warnings:
            toc_lines.append('- **[""Unexplained Errors"" #unexplained-errors]**')
            sections.append(('unexplained-errors', warnings))

        toc_lines.append('- **[""Info"" #info]**')
        sections.append(('info', self._get_general_info()))

        # Index of summary section.
        summary_index = len(sections)

        # Build a table containing summary functions of all other tables.
        # The actual section is added at position summary_index after creating
        # all other tables.
        if self.resolution in ['domain', 'combined']:
            summary = self._get_empty_table(title='summary')
            toc_lines.append('- **[""Summary"" #summary]**')

        for attribute in self.attributes:
            logging.info('Creating table(s) for %s' % attribute)
            tables = []
            if self.resolution in ['domain', 'combined']:
                if self.attribute_is_numeric(attribute):
                    domain_table = self._get_table(attribute)
                    tables.append(('', domain_table))
                    reports.extract_summary_rows(
                        domain_table, summary, link='#' + attribute)
                else:
                    tables.append((
                        '',
                        'Domain-wise reports only support numeric '
                        'attributes, but %s has type %s.' %
                        (attribute, self._all_attributes[attribute].__name__)))
            if self.resolution in ['problem', 'combined']:
                for domain in sorted(self.domains.keys()):
                    tables.append((domain, self._get_table(attribute, domain)))

            parts = []
            toc_line = []
            for (domain, table) in tables:
                if domain:
                    assert table
                    toc_line.append('[""%(domain)s"" #%(attribute)s-%(domain)s]' %
                                    locals())
                    parts.append('== %(domain)s ==[%(attribute)s-%(domain)s]\n'
                                 '%(table)s\n' % locals())
                else:
                    if table:
                        parts.append('%(table)s\n' % locals())
                    else:
                        parts.append('No task was found where all configurations '
                                     'have a value for "%s". Therefore no '
                                     'domain-wise table can be generated.\n' %
                                     attribute)

            toc_lines.append('- **[""%s"" #%s]**' % (attribute, attribute))
            toc_lines.append('  - ' + ' '.join(toc_line))
            sections.append((attribute, '\n'.join(parts)))

        # Add summary before main content. This is done after creating the main content
        # because the summary table is extracted from all other tables.
        if self.resolution in ['domain', 'combined']:
            sections.insert(summary_index, ('summary', summary))

        if self.resolution == 'domain':
            toc = '- ' + ' '.join('[""%s"" #%s]' % (attr, attr)
                                  for (attr, section) in sections)
        else:
            toc = '\n'.join(toc_lines)

        content = '\n'.join('= %s =[%s]\n\n%s' % (attr, attr, section)
                            for (attr, section) in sections)
        return '%s\n\n\n%s' % (toc, content)

    def _get_general_info(self):
        table = reports.Table(title='algorithm')
        for config, info in self.config_info.items():
            for attr in self.INFO_ATTRIBUTES:
                if info[attr]:
                    table.add_cell(config, attr, info[attr])
        table.set_column_order(self.INFO_ATTRIBUTES)
        return str(table)

    def _get_group_functions(self, attribute):
        """Decide on a list of group functions for this attribute."""
        return [(reports.function_name(f), f) for f in attribute.functions]

    def _add_table_info(self, attribute, func_name, table):
        """
        Add some information to the table for attributes where data is missing.
        """
        if not attribute.absolute:
            table.info.append('Only instances where all configurations have a '
                              'value for "%s" are considered.' % attribute)
            table.info.append('Each table entry gives the %s of "%s" for that '
                              'domain.' % (func_name, attribute))

        summary_names = [name.lower() for name, sum_func in table.summary_funcs.items()]
        if len(summary_names) == 1:
            table.info.append('The last row reports the %s across all domains.' %
                              summary_names[0])
        elif len(summary_names) > 1:
            table.info.append('The last rows report the %s across all domains.' %
                              ' and '.join(summary_names))

    def _get_suite_table(self, attribute):
        assert self.attribute_is_numeric(attribute), attribute
        table = self._get_empty_table(attribute)
        self._add_summary_functions(table, attribute)
        # The first group function is used for aggregation.
        func_name, func = self._get_group_functions(attribute)[0]
        num_probs = 0
        self._add_table_info(attribute, func_name, table)
        domain_config_values = defaultdict(list)
        for domain, problems in self.domains.items():
            for problem in problems:
                runs = self.problem_runs[(domain, problem)]
                if (not attribute.absolute and
                        any(run.get(attribute) is None for run in runs)):
                    continue
                num_probs += 1
                for run in runs:
                    value = run.get(attribute)
                    if value is not None:
                        domain_config_values[(domain, run['config'])].append(value)

        # If the attribute is absolute (e.g. coverage) we may have
        # added problems for which not all configs have a value. Therefore, we
        # can only print the number of instances (in brackets after the domain
        # name) if that number is the same for all configs. If not all configs
        # have values for the same number of problems, we write the full list of
        # different problem numbers.
        num_values_lists = defaultdict(list)
        for domain in self.domains:
            for config in self.configs:
                values = domain_config_values.get((domain, config), [])
                num_values_lists[domain].append(str(len(values)))
        for domain, num_values_list in num_values_lists.items():
            if len(set(num_values_list)) == 1:
                count = num_values_list[0]
            else:
                count = ','.join(num_values_list)
            link = None
            if self.resolution == 'combined':
                link = '#%s-%s' % (attribute, domain)
            formatter = reports.CellFormatter(link=link, count=count)
            table.cell_formatters[domain][table.header_column] = formatter

        for (domain, config), values in domain_config_values.items():
            table.add_cell(domain, config, func(values))

        table.num_values = num_probs
        return table

    def _get_domain_table(self, attribute, domain):
        table = self._get_empty_table(attribute)

        for config in self.configs:
            for run in self.domain_config_runs[domain, config]:
                table.add_cell(run['problem'], config, run.get(attribute))
        return table

    def _get_table(self, attribute, domain=None):
        if domain:
            return self._get_domain_table(attribute, domain)
        return self._get_suite_table(attribute)

    def _get_empty_table(self, attribute=None, title=None, columns=None):
        """Return an empty table."""
        if title is None:
            assert attribute is not None
            title = attribute
            if self.output_format == 'tex':
                title = title.capitalize().replace('_', ' ')
        if columns is None:
            columns = self.configs
        if attribute is not None and self.attribute_is_numeric(attribute):
            # Decide whether we want to highlight minima or maxima.
            min_wins = attribute.min_wins
            colored = self.colored and min_wins is not None
        else:
            # Do not highlight anything.
            min_wins = None
            colored = False
        table = reports.Table(title=title, min_wins=min_wins, colored=colored)
        table.set_column_order(columns)
        if self.resolution == 'combined':
            link = '#%s' % title
            formatter = reports.CellFormatter(link=link)
            table.cell_formatters[table.header_row][table.header_column] = formatter
        return table

    def _add_summary_functions(self, table, attribute):
        for funcname, func in self._get_group_functions(attribute):
            table.add_summary_function(funcname.capitalize(), func)
