====================================
DREAM: DiscRete sEmi-empiricAl Model
====================================

Authors:
Hao Fu <h.fu@soton.ac.uk>
Chris Marsden <c.marsden@soton.ac.uk>
Francesco Shankar <f.shankar@soton.ac.uk>


::: Prerequisites

* A python version 3.0 or above is recommended.
* Essential modules: 'numpy', 'scipy', 'matplotlib', 'colossus'
* Optional modules: 'joblib'

The folder "dream" contains the main code as well as
some functions that are called by the main code itself.

To compile the C-functions, type the following command:
    $ make clean
    $ make
or:
    $ make clean
    $ make _dream_

If 'make clean' returns the following messages
    $ rm dream/C_functions/*.so
    $ rm: dream/C_functions/*.so: No such file or directory
    $ make: *** [clean] Error 1
there is no need to worry, just proceed.


::: Run the code

Use the following command, for example, to call the code:
    $ python dream/DREAM.py run -o output -p explanatory.ini

Some commands are essential, without which the program will not run.
You need to precise the 'run' subcommand and the parameter file with '-p'.
The parameter file can have any extension. For information regarding the input
parameters see the detailed explanation in file "explanatory.ini".
You need also to precise the output folder name with '-o'. Non existing folder
will be created automatically.

The 'info' subcommand, meant to analyse outputs.
You need to specify the data folder as last argument of the command line.
Use the following command, for example, to have some info:
    $ python dream/DREAM.py info -p info_example.param output

Note this still under development.


::: Test run

To test the code is able to run correctly, launch a test run within a small
volume and mass range, say the parameters present in file "explanatory.ini".
You should get this print at terminal:

	Running DREAM (DiscRete sEmi-empiricAl Model)

	Cosmological model set to: planck18

	Reading parameters from file:
	path-to-DREAM/DREAM/explanatory.ini

	Generating catalog for a volume of (20.00 Mpc/h)^3

	Number of halos generated: 52

	Calculating accretion tracks...
	Accretion tracks calculated

	Generating halo IDs...
	Halo IDs generated

	Generating mergers...
	100%|██████████████████████████████████████████████████████████| 
	52/52 [00:00<00:00, 239.58it/s]

	Data stored in:
	path-to-DREAM/DREAM/outputs/data/

	Process completed correctly
	Time elapsed: 1 sec


::: Outputs

All output files will be store in a subfolder named "data",
inside the given folder, in txt format.
Existing files will be overwritten.

The output catalogue will be returned in two files:
1) output_parents.txt
2) output_mergers.txt

File 1) contains info on the parent haloes.
File 2) contains info on the subhaloes.

The output plots generated with subcommand 'info' will be stored in a subfolder
named "plots".