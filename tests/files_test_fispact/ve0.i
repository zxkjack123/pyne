<< CONTROL PHASE >>
<< enable JSON output >>
JSON 
<< overwrite existing output files of same name >>
CLOBBER 
<< perform collapse >>
GETXS 1 175
<< get decay data >>
GETDECAY 1
<< enable logging at level 1 >>
LOGLEVEL 1
<< approximate spectra when not available >>
SPEK 
<< end control >>
FISPACT 
* /home/gw/opt/pyne/tests/files_test_fispact/ve0

<< INITIALIZATION PHASE >>
<< output half life values >>
HALF 
<< output ingestion and inhalation values >>
HAZARDS 
<< set the target via FUEL >>
FUEL 2
H1 5.9740110000E+26
H2 6.8709030000E+22 
<< set the target density >>
DENSITY 1.0
<< set the threshold for atoms in the inventory >>
MIND 100000.0
<< output the initial inventory >>
ATOMS 

<< INVENTORY PHASE >>
<< irradiation schedule >>
FLUX 4060000000000.0
TIME 302400.0 SECS
ATOMS 
<< end of irradiation >>
FLUX 0.0
ZERO 
TIME 1.0 SECS
ATOMS 
<< end of cooling >>

<< end of input >>
END 
* end
