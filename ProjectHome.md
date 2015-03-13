# Introduction #

LDA-frames is an approach to identifying semantic frames from semantically unlabelled text corpora.
There are many frame formalisms but most of them suffer from the problem that all the frames must
be created manually and the set of semantic roles must be predefined. The LDA-Frames approach,
based on the Latent Dirichlet Allocation, avoids both these problems by employing statistics on
a syntactically annotated corpus. The only information that must be given is a number of semantic
frames and a number of semantic roles. This limitation, however, can be avoided
by automatic estimation of both these parameters. The model which is able to estimate the number of frames and roles automatically is called Non-parametric LDA-Frames.

More information and a demo application can be found on the <a href='http://nlp.fi.muni.cz/projekty/lda-frames/'>LDA-Frames homepage</a>.

# Dependencies #

  * Python >= 2.5
  * Numpy >= 1.6
  * NLTK >= 2.0
  * Scipy >= 0.9
  * Java 6 >= JDK 1.6 (Only when using the Stanford parser for generating the syntactic dependencies.)
  * g++ >= 4.6
  * Boost program options development package >= 1.46
  * GNU Scientific Library (GSL) development package >= 1.15
  * OpenMP (Only when using the multithreaded sampler)

# Building binaries from the source codes on Ubuntu 13.10. #

Install necessary dependencies.
```
$ sudo apt-get install g++ libboost-program-options-dev libgsl0-dev python-numpy python-nltk python-scipy
```

Move into the working directory of the sampler.
```
$ cd sampler
```

Build the binary file.
```
$ make
```

# Input and Output Data #

The input data for the sampler must have the form of a text file where each line corresponds to one lexical unit. Each line consists of tabulator-separated corpus realizations, where each realization is composed of space-separated positive integers that represent word identifiers for all slots. The number 0 is reserved for an empty slot. An example of the input file for 3 lexical units with two-slot realizations is shown below.

```
1 2	3 4	2 3
4 2	1 3	3 2	3 2
4 4	3 4	4 0
```

In order to be able to interpret the result, a dictionary that translates the line numbers into lexical units and word identifiers into words is required. The dictionary along with the sampler input file can be generated using the _generate_<u> </u>samplerinput.py_script. As its input, it takes a text file where each line
consists of a lexical unit and tabulator-separated realization from a corpus. The first line must start with ';' and must consist of tabulator-separated names of used grammatical relations. For example:_

```
;SUBJECT    OBJECT
eat 	people  	cake
eat 	man 		bread
cook    woman   	cake
```

The result produced by the sampler comprises two files – “frames.smpl”, which contains frame identifiers corresponding to the realization positions in the input file, and “roles.smpl”, which contains as many rows as and the number of frames, and as many columns as the number of slots. The numbers there represent semantic role identifiers.

# Example of the Usage #

Generate input data for the sampler along with a dictionary.
```
$ preprocessing/generate_samplerinput.py -l data/synthetic/synthetic.rel data/synthetic/
```

Run the Non-Parametric LDA-Frames sampler with 200 iterations.
```
$ sampler/ldaframes-sampler --seed=1 -I 100 data/synthetic/train.dat data/synthetic/
```

Evaluate the produced frames against the original data.
```
$ utils/evaluateRandom.py data/synthetic/synthetic.rel.chck data/synthetic/
```

Generate frame distributions for lexical units and word distributions for semantic roles.
```
$ preprocessing/generate_distributions.py data/synthetic/
```

Generate similatity matrix for all lexical units (based on the comparison of their frame distributions).
```
$ preprocessing/generate_similarities.py data/synthetic/
```

Show generated semantic frames for lexical unit "eat".
```
$ libldaf/showFrames.py -r 4 data/synthetic/ eat
```