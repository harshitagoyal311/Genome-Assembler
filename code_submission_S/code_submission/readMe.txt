ReadMe:

The file assembler V1 contains the code to generate contigs without using paired end data or the predictive model.

The file assembler V2 contains the code to generate contigs with the guidance of paired end data.

Files assemblerv1* and assemblerv2* have no external dependencies and can be directly compiled and executed as follow:

	g++ assembler<v1/v2>*.cpp <kmer-Size> <abundance>
	./a.out

File assembler V3 contains the code to generate contigs with the guidance of both paired end and machine learning model based combination of normalized scores. The code was tested on a OS X machine with NVIDIA tesla k80 GPU (on google cloud platform) and a lighter version on Satwik's machine (OS X NVIDIA 750M).


Before running the assembler file, we first need to get the trained model which is created inside the ipython-notebook file (which also includes the testing of other models such as logistic regression, random forest...)

To run the ipython notebook, follow these instructions:
	
	$ sudo easy_install pip				OSX
		or
	$ sudo apt-get install python-pip 	Linux

	$ pip install ipython[all]
	$ pip install tensorflow
	$ pip install sklearn
	$ ipython notebook

Go to the required directory and then open the notebook. The notebook takes the csv file generated from assembler V2 and creates the full branches (merged primary contig + branch string). The quality of each such branch is generated using CQAT and its formatted using clean_cqat_res.py file that can be found in the utils directory. The features for each contig is generated using the create_features.py file and/or the repeat_features.py file.

To clean CQAT results, run:
	
	$ python clean_cqat_res.py <contig.stats.tsv>

To generate features for each contig:

	$ python create_features.py <contig.fasta>

For our case the contig.fasta file contains each of the possible branches.


To install the dependencies for assembler V3 file, we first need to install bazel (google's api to build things related to Tensorflow) and protobuf (to load our trained graph in C++)

	$ brew install bazel

	$ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
	$ brew install protobuf

(for other OS, this manual can be referred: https://docs.bazel.build/versions/master/install.html)


We also need the tensorflow sources:

	$ mkdir /path/tensorflow
	$ cd /path/tensorflow
	$ git clone https://github.com/tensorflow/tensorflow.git

	$ cd /path/tensorflow
	$ ./configure

Once the dependencies have been install, the code can built and run by using bazel and the build file:

	$ bazel run -c opt //tensorflow/cc/models:model