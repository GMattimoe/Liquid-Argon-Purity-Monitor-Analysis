Liquid argon peak amplitude project

These programs are written in c++ and ROOT.
It is necessary to have c++ and ROOT installed and the FFTW Wrapper in the 
path. Voltage values were recorded from the cathode and anode of a liquid argon
purity monitor with their corresponding time. These voltage measurements were 
recorded in sets of 1000 at a single time. These sets were taken every 30 
minutes over a period of days. This raw data was smoothed and averaged and 
the the cathode averages were given unique names ending in 
".ch4.traces_averages.root" with times, ".ch4.traces.times"

newCathodeTimePlot.cpp
This program loops over all the ".ch4.traces_averages.root" files in the 
directory and extracts the peak amplitude, mapping it to the time in 
the corresponding ".ch4.times" file. The program writes the TGraph of peak 
amplitude agains time to a .root file

fittingTime.cpp
This progam fits the peak amplitude graph with a number of functions and
prints the graphs.