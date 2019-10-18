# panopticon

![Utopian data structure](https://upload.wikimedia.org/wikipedia/en/e/e1/Panopticon_Willey_Reveley_1791.png )



Panopticon is a set of tools named after Jeremy Bentham's model for a prison, whereby multiple cells (!) could be observed with ease by a single individual.

You can install the panopticon on the command line by cloning this git repository :
```
    git clone https://smarkson@bitbucket.org/carterlab/panopticon.git
    cd panopticon
    pip install --editable .
    pyensembl install --release 75 --species homo_sapiens
```

Now try exploring the exciting panopticon programming with

```
    panopticon --help
```

If you would like autocompletion of panopticon commands, put the following line in your `.bashrc`
```
eval "$(_PANOPTICON_COMPLETE=source panopticon)"
```

# To do
- Clean up docstrings and make sure that they are in numpydoc format
- add in sphinx documentation
- add in tests, including test data
- add in test of signature capability
- add in difference of correlations test
