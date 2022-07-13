# pywrap for axion-fermion stars

Example:
```bash
python3 wrapper.py par/axion_fa_exp_-1d7.par fa-1.7 &
```
Compile and run Shoot + Solve using the input written in the parfile `axion_fa_exp_-1d7.par`.
Create the file `log_wrappy_${localtime}.py` where info on the run are printed. The output is
saved in the directory `./out/fa-1.7`. If the directory already exists, the program stops.
If axion-star (i.e. `rho0c=0`), then create also a M(phic) plot: `./out/fa-1.7/existence_plot.png`.

For compiling the level-curver script, go in `mass_curves/`and then execute:
```bash
ln -s Examples/AnalyticalExample.cpp
g++ -o example AnalyticalExample.cpp LevelCurves.cpp
./example
```

This will produce a file `out_example.txt` that contains a level curve.
To plot the result, you can open `gnuplot` and use
```bash
plot 'out_example.txt' u 2:1 i 1 w l
```

You can check the result against the output of the `Matlab` script `LevelCurves.m`:
```bash
x_init          = -0.6;
y_init          = -0.30103;
start_direction = 2;
ds              = 0.01;
maxpoints       = 100;
onlySquares     = 0;
saveGIF         = 0;
LevelCurves(x_init, y_init, start_direction, ds, maxpoints, onlySquares, saveGIF)
```


