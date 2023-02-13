#!/usr/bin/env python3

# This script wrangles a csv-table of author names into a correctly formatted author list for journal publication.
# Written by Johannes Englisch.

# python3 wrangle_author_list.py author_list.csv > author_list.html

from collections import OrderedDict
import re
import sys

from csvw import dsv


AFFILIATION_COLS = [
   'First affiliation',
   'Second affiliation',
   'Third affiliation',
   'fourth affiliation',
   'fifth affiliation',
]

NOTE_COLUMN = 'Author notes description'
NOTE_SYMBOLS = ['*', '†', '‡', '§', '¶', '#']

HTMLDOC = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
</head>
<body>
{}
</body>
</html>"""


def note_symbol(note_index):
   return (
       NOTE_SYMBOLS[note_index % len(NOTE_SYMBOLS)]
       * (note_index // len(NOTE_SYMBOLS) + 1))


def name_to_initial(name):
   return '-'.join('{}.'.format(n[0]) for n in name.split('-'))


def initials(first_names):
   return ' '.join(map(name_to_initial, first_names.split()))


def format_author(person, people_affiliations, people_notes):
   affiliations = ','.join(
       map(str, people_affiliations.get(person['index'], ())))
   notes = ''.join(
       map(note_symbol, people_notes.get(person['index'], ())))

   return '{first_name} {last_name}{affils}{notes}'.format(
       first_name=person['first name'],
       last_name=person['last name'],
       affils='<sup>{}</sup>'.format(affiliations) if affiliations else '',
       notes=notes)


def main():
   try:
       file_name = sys.argv[1]
   except IndexError:
       print('usage:', sys.argv[0], 'filename', file=sys.stderr)
       sys.exit(1)

   # read data

   rows = dsv.reader(file_name, encoding='utf-8', dicts=True)
   rows = [
       OrderedDict(
           (k.strip(), v.strip())
           for k ,v in row.items()
           if k and v and k.strip() and v.strip())
       for row in rows]
   rows = list(filter(None, rows))

   # collect affiliations

   affiliation_numbers = OrderedDict()
   note_indices = OrderedDict()
   # predefined by Science
   note_indices['Corresponding author'] = 0

   people_affiliations = OrderedDict()
   people_notes = OrderedDict()

   for person in rows:
       person_index = person['index']

       for col in AFFILIATION_COLS:
           affiliation = person.get(col)
           if affiliation and affiliation not in affiliation_numbers:
               affiliation_numbers[affiliation] = len(affiliation_numbers) + 1
       person_affiliations = [
           affiliation_numbers[person[col]]
           for col in AFFILIATION_COLS
           if col in person]
       if person_affiliations:
           people_affiliations[person_index] = person_affiliations

       notes = [s.strip() for s in person.get(NOTE_COLUMN, '').split(';')]
       if notes:
           for note in notes:
               if note and note not in note_indices:
                   note_indices[note] = len(note_indices)
           people_notes[person_index] = [
               note_indices[note]
               for note in notes
               if note]

   # output

   print(HTMLDOC.format('\n'.join((
       '<p>',
       ',\n'.join(
           format_author(person, people_affiliations, people_notes)
           for person in rows),
       '</p>',

       '<p><strong>Affiliations:</strong></p>',

       '\n'.join(
           '<p><sup>{}</sup> {}</p>'.format(number, affiliation)
           for affiliation, number in affiliation_numbers.items()),

       '<p><strong>Notes:</strong></p>',

       '\n'.join(
           '<p>{}{}</p>'.format(note_symbol(i), n)
           for n, i in note_indices.items())))))


if __name__ == '__main__':
   main()
