============
Installation
============

You can install the panopticon on the command line by cloning this git
repository. Because this package is in a rather alpha mode, I recommend
installing with the ``--editable`` flag.

::

        git clone https://github.com/scyrusm/panopticon.git
        cd panopticon
        pip install --editable .

For some functions you will also need to run

::

        pyensembl install --release 75 --species homo_sapiens

Now try exploring the exciting panopticon programming with

::

        panopticon --help

If you would like autocompletion of panopticon commands (though I don't
do this myself), put the following line in your ``.bashrc``

::

    eval "$(_PANOPTICON_COMPLETE=source panopticon)"
