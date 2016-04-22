# Superpixel Hierarchical Clustering algorithm (SPHC) For Image Segmentation

The idea to create the algorithm came from difficulties in locating online resources describing image segmentation algorithms that use superpixels as a starting point. The conception of the algorithm further came from the observation that neighboring superpixels often have similarities in color and object boundaries are defined by color differences/similarities.

The algorithm takes two main inputs: a RGB pixel grid to represent an image and a grid of segments from the sklearn SLIC superpixel-creating algorithm. 

After segment assignment happens, the superpixel hierarchical clustering takes place as follows
1. Loop through these four steps while the smallest cluster euclidean distance is below a specified threshold:
1a. Find neighboring segment pairs for each of the 1 through K superpixels.
1b. Get average R, G, and B values for each segment
1c. For each pair of neighboring segments, calculate euclidean distance using average R, G, and B values.
1d. Merge the 2 segments with the shortest RGB euclidean distance.
2. Output image.

##References
1. Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aure-lien Lucchi, Pascal Fua, and Sabine Susstrunk, SLIC Superpixels, EPFL Technical Report 149300, June 2010.
2. Adrian Rosebrock, A SLIC Superpixel Tutorial using Python, http://www.pyimagesearch.com/2014/07/28/a-slic-superpixel-tutorial-using-python, July 28, 2014.
