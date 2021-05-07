ADS to BibDesk
==============

This is the command line edition of ADS to BibDesk, a tool for retrieving the bibtex, abstract and PDF of an astronomical journal article published on `ADS <http://adsabs.harvard.edu>`_ or `arXiv.org <http://arxiv.org/archive/astro-ph>`_ and add it to your `BibDesk <http://bibdesk.sourceforge.net/>`_ database.

ADS to BibDesk is a tool for retrieving the bibtex, abstract and PDF of an astronomical journal article published on `ADS <http://adsabs.harvard.edu>`_ or `arXiv.org <http://arxiv.org/archive/astro-ph>`_ and adding it to your `BibDesk <http://bibdesk.sourceforge.net/>`_ database.

ADS to BibDesk comes in two flavours: an Automator Service that you can use to grab papers in any app (e.g., in Safari, or Mail), or a command line app.

*Developers:* please read the `CONTRIBUTING <https://github.com/jonathansick/ads_bibdesk/blob/master/CONTRIBUTING.md>`_ document for details on how to build the ADS to BibDesk CLI/Service from source, make changes, and submit pull requests.

Command Line Quickstart
-----------------------

ADS to BibDesk can also be run directly from the command line.
The command line script can be installed via::

    python setup.py install

You may need to run the last command with `sudo`.

Once `adsbibdesk` is installed, you can call it with the same types of article tokens you can launch the Service with, e.g.,::

    adsbibdesk 1998ApJ...500..525S

A full summary of `adsbibdesk` commands is available via::

    adsbibdesk --help

Summary of article tokens
-------------------------

* The URL of an ADS or arXiv article page,
* The ADS bibcode of an article (e.g. `1998ApJ...500..525S`),
* The arXiv identifier of an article (e.g. `0911.4956`), or
* An article DOI.

Other Modes
-----------

Besides the primary mode (adding a single paper to BibDesk, ADS to BibDesk has three other modes: previewing papers, updated preprints, and ingesting PDF archives into BibDesk.

Previewing Papers
~~~~~~~~~~~~~~~~~

Use the `-o` switch to simply download and view the PDF of an article without adding it to BibDesk. E.g.,::

    adsbibdesk -o 1998ApJ...500..525S

Updating Preprints
~~~~~~~~~~~~~~~~~~

Run ADS to BibDesk with the `-u` switch to find and update all astro-ph preprints in your BibDesk bibliography::

    adsbibdesk -u

To restrict this update to a date range, you can use the `--from_date` (`-f`) and `--to_date` (`-t`) flags with dates in `MM/YY` format. For example, to update preprints published in 2012, run::

    adsbibdesk -u --from_date=01/12 --to_date=12/12

Note that this operation can take some time since we throttle requests to ADS to be a nicer robot.

PDF Ingest Mode
~~~~~~~~~~~~~~~

With the command-line ADS to BibDesk, you can ingest a folder of PDFs that originated from ADS into BibDesk.
This is great for users who have amassed a literature folder, but are just starting to use BibDesk.
This will get you started quickly.

You need the program `pdf2json <http://code.google.com/p/pdf2json/>`_ to use
this script. The easiest way to get pdf2json and its dependencies is through
`Homebrew <http://mxcl.github.com/homebrew/>`_, the Mac package manager.
Once homebrew is setup, simply run `brew install pdf2json`.

To run this workflow,::

    adsbibdesk -p my_pdf_dir/

where `my_pdf_dir/` is a directory containing PDFs that you want to ingest.

Note that this workflow relies on a DOI existing in the PDF.
As such, it will not identify astro-ph pre-prints, or published papers older than a few years.
Typically the DOI is published on the first page of modern papers.
This method was inspired by a script by `Dr Lucy Lim <http://www.mit.edu/people/lucylim/BibDesk.html>`_.

License
-------

Copyright 2014 Jonathan Sick, Rui Pereira and Dan-Foreman Mackey

ADS to BibDesk is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ADS to BibDesk is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ADS to BibDesk.  If not, see <http://www.gnu.org/licenses/>.
