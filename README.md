# idl_mommaps
IDL package for generating moment maps from radio data cubes

Requires:


[IDL Astronomy Library:](http://idlastro.gsfc.nasa.gov)
        <https://github.com/wlandsman/IDLAstro>

[Coyote Graphics Plotting Library:](http://www.idlcoyote.com/documents/programs.php)
        <https://github.com/idl-coyote/coyote>


The main program to run is [*makemom.pro*](https://github.com/tonywong94/idl_mommaps/blob/master/makemom.pro).  Input parameters are summarized in the header of that program.

As of 2020, this package is no longer being developed - please switch to the [python version](https://github.com/tonywong94/maskmoment) instead.

Contributions from Tony Wong, Rui Xue, Annie Hughes, and Erik Rosolowsky.

* v1 - 27may2015 - initial release, still cleaning up documentation.
* v2 - 04jun2015 - include plotting routines
* v3 - 02feb2016 - new approach to estimating uncertainties in mom0, flux
