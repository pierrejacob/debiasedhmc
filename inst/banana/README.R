This folder has scripts that implements the banana shaped distribution.
On this example we implement the contractive coupling presented in
Bou-Rabee, Eberle and Zimmer 2018 https://arxiv.org/abs/1805.00452
and observe that it leads to shorter meeting times.

The script banana.synchronous.run.R obtains meeting times
by using a synchronous coupling, i.e. common random numbers
for the initial momentum variables.

The script banana.contractive.run.R obtains meeting times
by using a contractive coupling, i.e. using a mix of reflection
coupling and maximal coupling for the initial momentum variables.

Lastly, banana.plots.R produces histograms of the meeting times.


