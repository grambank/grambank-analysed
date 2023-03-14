#!/usr/bin/env python3

# This script wrangles a csv-table of author names into a correctly formatted
# bibtex entry.
# Written by Johannes Englisch.

# python3 make_bibtex.py author_list.csv > grambank-entry.bib

from collections import OrderedDict
import csv
from itertools import islice
import unicodedata
import sys


# # TODO: add volume and number once they're known
# volume = {{{volume}}},
# number = {{{number}}},
BIBTEX_TEMPLATE = """@{bibtype}{{{bibkey},
  author = {{{author}}},
  year = {{{year}}},
  title = {{{title}}},
  journal = {{{journal}}},
}}
"""


NON_ASCII = set()

REPLACEMENTS = [
    ('ß', r'{\ss}'),
    ('á', r"{\'a}"),
    ('â', r'{\^o}'),
    ('å', r'{\aa}'),
    ('é', r"{\'e}"),
    ('ö', r'{\"o}'),
    ('ü', r'{\"u}'),
]


def latexify_name(name):
    name = name.strip()
    name = unicodedata.normalize('NFKC', name)
    for old, new in REPLACEMENTS:
        name = name.replace(old, new)
    return name


def format_author(row):
    first_name = latexify_name(row['first name'])
    last_name = latexify_name(row['last name'])
    if '葉婧婷' in last_name:
        return r'Ye \begin{CJK*}{UTF8}{bsmi}葉婧婷\end{CJK*}, Jingting'
    else:
        NON_ASCII.update(char for char in first_name if ord(char) > 127)
        NON_ASCII.update(char for char in last_name if ord(char) > 127)
        return '{}, {}'.format(last_name, first_name)


def main():
    try:
        file_name = sys.argv[1]
    except IndexError:
        print('usage:', sys.argv[0], 'filename', file=sys.stderr)
        sys.exit(1)

    # read data

    with open(file_name, encoding='utf-8') as f:
        rows = list(csv.reader(f))

    header = [cell.lstrip('\ufeff').strip() for cell in rows[0]]
    # strip trailing empty cols
    while header and not header[-1].strip():
        header.pop()

    rows = [
        OrderedDict(
            (k, v.strip())
            for k, v in zip(header, row)
            if v and v.strip())
        for row in islice(rows, 1, len(rows))]
    rows = [row for row in rows if 'index' in row]

    # format authors

    authors = ' and '.join(map(format_author, rows))
    if NON_ASCII:
        print(
            'Non-ASCII characters:',
            ', '.join(sorted(NON_ASCII)),
            file=sys.stderr)

    # output

    # volume='???',
    # number='???',
    print(BIBTEX_TEMPLATE.format(
        bibtype='article',
        bibkey='Grambank2023',
        author=authors,
        year='in press',
        title='Grambank reveals the importance of genealogical constraints on linguistic diversity and highlights the impact of language loss',
        journal='Science Advances',
    ))


if __name__ == '__main__':
    main()
