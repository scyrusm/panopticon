============
Installation
============
Panopticon can be installed from pypi using:
::

       pip install panopticon-single-cell

::


Or you can install panopticon on the command line by cloning this git
repository. Because this package is in active development, I recommend
installing with the ``--editable`` flag.

::

        git clone https://github.com/scyrusm/panopticon.git
        cd panopticon
        pip install --editable .

For some functions you will also need to run (release subject to change)

::

        pyensembl install --release 75 --species homo_sapiens

Now try exploring the exciting panopticon programming with

::

        panopticon --help

If you would like autocompletion of panopticon commands (though I don't
do this myself), put the following line in your ``.bashrc``

::

    eval "$(_PANOPTICON_COMPLETE=source panopticon)"

File locking
============

I strongly recommend using file locking, to ensure that you do not accidentally open a .loom file outside of a 'with' context in two place concurrently.  If you do this, your .loom file may become corrupted.  

To ensure that file locking is used, set the following environmental variable, e.g. in your ``.bashrc``

::

    export HDF5_USE_FILE_LOCKING=TRUE
