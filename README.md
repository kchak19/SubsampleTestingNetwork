README File



Numerical Studies:

***** Note the following ******
	Parallel computation is necessary for the subsampling-based methods. Using a large number of cores can be significantly beneficial in terms of time saving. However, each core needs to access enough memory (RAM) depending on subsample sizes.
	In all numerical studies we demonstrate our advantage using time elapsed. However, it can vary based on different computing machines.
	Since the testing method without subsampling can take a huge amount of time, we have been conservative in choosing several aspects such as the number of bootstrap resamples, the number of iterations and so on. The user can attempt to run the testing method without subsampling as long as they understand the time concerns.
	Graph sizes and subsampling parameters for the simulation experiments are kept at low values in the codes than in the paper. The user can change them as required.


We keep d=2 (the rank of RDPG model) for all the simulations, as well as real data examples. 

Simulation Experiment for Estimation: (FinEst.R)

Specify the number of nodes n and the vector of subsampling parameter values ms and ks. The number of replications to account for randomness, as described in the paper, is denoted by repl. This experiment will give outputs of mean relative norm difference errors of the without subsampling and with subsampling estimation methods, as well as the variations and average times elapsed. This is how Tables 1 and 2 are generated in the paper.

Simulation Experiment for Testing: (FinTest.R)

Specify n,k,m, and the number of bootstrap resamples bs used for testing. We calculate the level and power of the tests without subsampling, with subsampling and multiple testing. The eps variable denotes the ϵ that has been used in the paper to indicate the difference between the null and the alternative. Changing the ϵ to different values gives us Table 3 presented in the paper. Note that the level corresponds to power at ϵ=0. Runtimes for these methods will also be generated as output.


Real Data: FriendFeed Data (FriendFeed.R)

We consider undirected, unweighted networks with no self-loops consisting of three layers: likes,comments and follows. These networks are compared pairwise using bootstrap testing with subsampling and multiple testing. For convenience, only one pairwise comparison among the three is active, and others are commented out. Only one choice of subsampling parameters is used, and it can be replaced by other parameter choices. The testing without subsampling can be easily included using the appropriate functions (setting  method="ASE" in testNetwork). This is how Table 4 is generated in the paper.

Real Data: Sanremo Data (Sanremo.R)

Based on the Sanremo Network data, choose only nodes that have high weight index (highWInd) and consider undirected, unweighted networks with no self-loops consisting of three layers: retweets,mentions and replies. These networks are compared using bootstrap testing with subsampling and multiple testing. For convenience, only one pairwise comparison among the three is active, and others are commented out. Only one choice of subsampling parameters is used, and it can be replaced by other parameter choices. The testing without subsampling can be easily included using the appropriate functions (setting method="ASE" in testNetwork). This is how Table 4 is generated in the paper.


In the following, we list all the functions that we have used for the numerical studies.

Functions: (Functions.R)


	graph(P)
	Objective: Generates a graph from a probability matrix.
	Input: 
	P is the probability matrix from which the graph is generated.
	Output: 
	An adjacency matrix of the same size as P.

	ASE(A,dim)
	Objective: Computes the Adjacency Spectral Embedding (ASE) on a fixed dimension.
	Input: 
	A is an adjacency matrix.
	dim is the fixed dimension for ASE.
	Output:
	The ASE of matrix A at dimension dim.

	procr(X,Y)
	Objective: Finds the Procrustes rotation between X and Y.
	Input:
	Two matrices X and Y of the same dimension.
	Output:
	The minimum Frobenius norm distance between X and Y.
	The orthogonal rotation between X and Y.

	ProbM(X)
	Objective: Computes the probability matrix from a latent matrix.
	Input:
	Latent matrix X.
	Output:
	Probability matrix.

	estSS(A,m,k,d,overlap)
	Objective: Finds the ASE based estimated latent matrix on a fixed dimension with subsampling.
	Input:
	A is an adjacency matrix.
	m is the overlap region size.
	k is the number of subsamples.
	d is the fixed dimension for ASE.
	overlap is the method for choosing the overlap region. Default is "dense". Another choice is "random".
	Output:
	The ASE based estimated latent matrix of matrix A at dimension d.
	Note: This function makes use of parallel computation. It is recommended to define appropriate parallel clusters with enough cores to get better results.

	EST(A,m,k,d,method)
	Objective: Finds the estimated latent matrix of A based on different methods.
	Input:
	A is an adjacency matrix.
	m is the overlap region size.
	k is the number of subsamples.
	d is the fixed dimension for ASE.
	method is the method for performing estimation. Default is "SS" (subsampling). Another choice is "ASE".
	Output:
	The estimated latent matrix of matrix A at dimension d.
	Note: This is a function combining ASE and estSS, created for convenience in later simulation experiments. The arguments m,k are unnecessary when using method="ASE", and they can be specified arbitrarily.

	testNetwork(A,B,m,k,d,bs,method)
	Objective: Performs bootstrap based testing between two adjacency matrices of two networks.
	Input:
	A is an adjacency matrix.
	B is another adjacency matrix.
	m is the overlap region size.
	k is the number of subsamples.
	d is the fixed dimension for ASE.
	bs is the number of Bootstrap resamples.
	method is the method for performing estimation. Default is "SS" (subsampling). Another choice is "ASE".
	Output: A list containing:
	pval is the p-value for the test.
	time is the time taken for the method in performing the bootstrap test.
	Note: To perform the full network bootstrap testing method (ASE), we do not need the values m,k, and they can be specified arbitrarily.

	multTestNetwork(A,B,m,k,d,bs)
	Objective: Performs multiple testing between two adjacency matrices of two networks.
	Input:
	A is an adjacency matrix.
	B is another adjacency matrix.
	m is the overlap region size.
	k is the number of subsamples.
	d is the fixed dimension for ASE.
	bs is the number of Bootstrap resamples.
	Output: A list containing:
	pval is the p-value for the test.
	time is the time taken for the method in performing the multiple testing.
	Note: For each of the tests on the subsamples in the multiple testing procedure, ASE based full network estimation is used.

![image](https://github.com/user-attachments/assets/6285793c-59d2-4e62-8e7a-c733b86f6914)
