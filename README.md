The aim of libgpr is to make using Genetic Programming easy to include within any C/C++ application.  Genetic programming (GP) is a powerful technique, inspired by the process of natural selection, which can be utilized to automatically discover programs which produce a desired input to output transformation.  The problem can be stated simply as:

"**Given some behavior find a program which produces this behavior, or a close approximation of it**"

Both classical tree-based and Cartesian forms of Genetic Programming are supported, including self-modifying CGP.  It's also possible to export the results in the form of C programs suitable for running on Arduino micro-controllers. Provided that the simulation system used to evaluate fitness is somewhat similar to a real hardware implementation then you should be able to upload the program to an Arduino, plug in some inputs and outputs and obtain similar types of behavior.


Installation
============

To build from source:

```bash
    sudo apt-get install gnuplot libz-dev
    make
    sudo make install
```

This creates the library and installs it into /usr/local.  If you need to run the unit tests after installing as above:

```bash
    make tests
    ./libgpr_tests
```

Or to check for any memory leaks:

```bash
    valgrind --leak-check=full ./libgpr_tests
```

There are also some example programs within the libtest and libtest_cartesian directories, which can be built by running:

```bash
    make ltest
    ./libgpr
```

or

```bash
    make ltestc
    ./libgpr
```

Usage
=====

For more detailed information on usage see http://robotics.uk.to/doku.php?id=libgpr or view the manpage.

References
==========

For more information on the Genetic Programming method, see:

    * Genetic Programming: On the programming of computers by means of natural selection, John R. Koza, MIT Press

    * https://en.wikipedia.org/wiki/Genetic_programming

    * http://www.cartesiangp.co.uk

    * https://en.wikipedia.org/wiki/Arduino

    * Self-Modifying Cartesian Genetic Programming, Simon Harding

    * Automatic Algorithm Invention with a GPU  https://www.youtube.com/watch?v=xQDazGrKsuM
