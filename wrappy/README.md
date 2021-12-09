# pywrap for axion-fermion stars

Example:
```bash
python3 wrapper.py par/axion_fa_exp_-1d7.par fa-1.7 &
```
Compile and run Shoot + Solve using the input written in the parfile `axion_fa_exp_-1d7.par`.
Create the file `log_wrappy_${localtime}.py` where info on the run are printed. The output is
saved in the directory `./out/fa-1.7`. If the directory already exists, the program stops.
If axion-star (i.e. `rho0c=0`), then create also a M(phic) plot: `./out/fa-1.7/existence_plot.png`.

To see if the code is running type
	top
To stop he code, type
	kill NUMBER_OF_PROCESS_FROM_TOP_COMAND
