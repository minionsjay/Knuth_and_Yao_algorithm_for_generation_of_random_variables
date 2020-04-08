# Knuth_and_Yao_algorithm_for_generation_of_random_variables

This a C++ implementation of an algorithm from Knuth and Yao published in 1976 originally. The algorithm uses a source of unbiased i.i.d. bits to generate a random variable according to a probability distribution stored as a probability vector. It uses the binary expansions of the probability values to create a tree (DDG tree) to generate a random variable. If H denotes the entropy of the given distribution, then the algorith is optimal from an information theoritical point of view as it uses an expected number of bits between H and H + 2.

The code only uses the standard library. Some binary files describing a few examples from the PDF file can be downloaded from here too.

Please refer to the PDF "ky_algo.pdf" for more information and references.

The code uses NTL ( https://www.shoup.net/ntl/ ) to compute p-bit accurate probability values.
