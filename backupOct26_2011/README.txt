How to add a new algorithm and plot results

1) Edit the script processdat.m and add your algorithm at the end of the list

allalgs(48).algo = 'matlab_function';   % name of Matlab function
allalgs(48).speed = 0;                  % speed of algorithm (simply enter 0)
allalgs(48).name = 'This my ICA';       % name of Algorithm (for figures etc...)
allalgs(48).options = { };              % options to give to the Matlab function (in addition to the data)

You may add several ICA algorithms by incrementing the index to 49, 50, etc...

2) Your algorithm has index 48, to run it on a given dataset use

DATASET=1; ALGONUM=48; processdat;
DATASET=2; ALGONUM=48; processdat;
DATASET=3; ALGONUM=48; processdat;
DATASET=4; ALGONUM=48; processdat;
DATASET=5; ALGONUM=48; processdat;
DATASET=6; ALGONUM=48; processdat;
DATASET=7; ALGONUM=48; processdat;
DATASET=8; ALGONUM=48; processdat;
DATASET=9; ALGONUM=48; processdat;
DATASET=10; ALGONUM=48; processdat;
DATASET=11; ALGONUM=48; processdat;
DATASET=12; ALGONUM=48; processdat;
DATASET=13; ALGONUM=48; processdat;
DATASET=14; ALGONUM=48; processdat;

3) Compute MIR (mutual information reduction) for all datasets

In the function mutualinfoalgo.m

Add the name of your algorithm (MUST be the same name as in allalgs(48).name)
to the list of algorithm at the beginning of the function.
The function will rescan all algorithm and datasets. It may take several
hours.

4) Add the name of your algorithm to the beginning of the "plotresults.m"
script. Because some ICA algorithm are computed and not plotted
you must add the name of your algorithm to the list.

************************************
Arnaud Delorme - September 1st, 2010