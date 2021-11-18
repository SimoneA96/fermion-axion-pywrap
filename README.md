# pywrap for axion-fermion stars

Example:
```bash
python3 wrapper.py par/axion_fa_exp_-1d7.par fa-1.7 &
```
Compile and run Shoot + Solve using the input written in the parfile `axion_fa_exp_-1d7.par`.
Create the file `log_wrappy.py` where info on the run are printed. The output is
saved in the directory `./out/fa-1.7`. If axion-star (i.e. `rho0c=0`), then
create also a M(phic) plot: `./out/fa-1.7/existence_plot.png`
