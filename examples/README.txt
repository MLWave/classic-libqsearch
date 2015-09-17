The best examples for this package are in the libcomplearn-dev package:
/usr/share/doc/libcomplearn-dev/examples

To use the directories, first make a distance matrix, for example:
To make an output CompLearn Binary distance matrix, use the -b option.
ncd -d 10-mammals 10-mammals -b
then, use maketree:
maketree distmatrix.clb

