This repo contains software to accompany the paper, "Empirically Measuring Online Social Influence". The contents is as follows:

	- *00_raw*: This folder contains the tweets ids required to rehydrate the dataset, as well as the final empirical influence scores.
	- *01_simulation*: This folder contains a script responsible for simulating the measurement framework. It uses the R `future` package to run simulations across an ad hoc cluster of nodes.
	- *02_mturk*: This folder contains the templates for various experimental contexts and a script that runs the empirical measurement framework. Note that MTURK environment variable setup is required for this to work.
