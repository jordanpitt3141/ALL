set table "Wavetank.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set samples 25; set dummy x; plot [x=1:10] plot 'file.dat' smooth cspline; ;
