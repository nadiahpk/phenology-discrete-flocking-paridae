# Brief Description

Octave code to find the evolutionarily singular breeding-time strategy for a strict, discrete-flocking resident bird species (e.g. marsh tit) in a phenologically heterogenous landscape.

# Reference

Kristensen, N.P., Johansson, J., Smith, H.G., Jonzen, N. Paper in preparation.

# Quick Start

If you have git:

    $ git clone https://github.com/nadiahpk/phenology-discrete-flocking-paridae
    $ cd phenology-two-trait-migratory-bird

If you don't have git, just download and unpack the latest zip file:

[https://github.com/nadiahpk/phenology-discrete-flocking-paridae/archive/master.zip](https://github.com/nadiahpk/phenology-discrete-flocking-paridae/archive/master.zip)

In Octave,

    > marsh2

will load the parameter dictionary ```p``` into the workspace. You can then define a good initial guess for the system's evolutionarily singular state. Our first guess for the hatching dates in each habitat type is

    > x0=[124   164]

and our first guess for the population size is

    > n0 = [49      49]

The solver can then be run,

    > x=calc_x(p,x0,n0)

    x = 124.33   165.27

    > n=calc_n(p,x,n0)

    n = 49.217   49.252

and this is our singular strategy hatching dates and population size respectively.

## Issues

The solver needs a fairly good initial guess in order to succeed at finding the singular strategy. Some helpful tips:

* Prior residence means that the hatching date will always be earlier than synchrony, so make sure that ```x0``` is lower than ```p.x_opt```. The solver can get "lost" if you start at the wrong end.
* The biologically reasonable system will be close to saturated, so you can often use ```K``` as a good initial guess for ```n```.
* Do you already know the singular strategy nearby in the parameter space? Perhaps start there and work your way slowly to the area of interest. If you are sweeping the parameter space, it's a good idea to go through incrementally, using the mismatch of the singular strategy at the previous step as the initial guess for the next step.

# Taking it further

(to be written)

* Evolutionary stability analysis
* Systems of three or more habitat types
