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

## Find an evolutionarily singular strategy

In Octave,

```> default_parameter_values```

will load the parameter dictionary p into the workspace. 

The default parameters have both habitats with the same
plant phenology

```
> p.x_opt
171   171
```

Let's change it so that the first habitat has a phenology
five days earlier than the second habitat

```
> p.x_opt(1) = 171-5;
> p.x_opt
ans =

   166   171 
```

You can then define a good initial guess for the system's
evolutionarily singular state. Let's say our first guess for the
hatching dates in each habitat type is about 8 days earlier
than the optimal hatching date

```
> x0 = [159 164];
```

and our first guess for the population size we'll set equal
to the number of available territories in each habitat type

```
> n0 = p.K
n0 =

  50   50 
```

The solver can then be run,

```
> x=calc_x(p,x0,n0)

x = 158.81   163.11

> n=calc_n(p,x,n0)

n = 46.764   46.692
```

and this is our singular strategy hatching dates and
population size respectively.

## Issues

The solver needs a fairly good initial guess in order to succeed at finding the singular strategy. Some helpful tips:

* Prior residence means that the hatching date will always be earlier than synchrony, so make sure that ```x0``` is lower than ```p.x_opt```. The solver can get "lost" if you start at the wrong end.
* The biologically reasonable system will be close to saturated, so you can often use ```K``` as a good initial guess for ```n```.
* Do you already know the singular strategy nearby in the parameter space? Perhaps start there and work your way slowly to the area of interest. If you are sweeping the parameter space, it's a good idea to go through incrementally, using the mismatch of the singular strategy at the previous step as the initial guess for the next step.

# Taking it further

## Check the evolutionary stability of a singular strategy

Let's check the evolutionary stability of the singular
strategy above. In the same workspace

```
> [eigH,eigJ,eigJs,Hess,Jac] = check_stab(p,x,n);
> eigH
eigH = -0.0012491
> eigJs
eigJs = -0.0016523
```

The dominant eigenvalue of the Hessian matrix ```eigH``` is
negative, so this singular strategy is a fitness maximum.
The dominant eigenvalue of Leimar's (2009) Jacobian
matrix Js is negative, so the singular strategy is also
an evolutionary attractor. Therefore we have found a singular
strategy that is a continuously stable strategy (CSS).

## Reproduce Figure 2a

First we'll read in the default parameter values supplied 

```
> default_parameter_values
```

We want to vary the landscape-composition parameter ```K```.
It is stored in the parameter values dictionary ```p``` as a
vector of two values. We'd like to vary it between 20 and 80
percent early habitat, so we create a matrix with each row
corresponding to a point at which we'd like to solve the
system, and the columns being early- and late-K
```
> parV = [linspace(20,80,13)',linspace(80,20,13)'];
```

We need to start the system off with an initial guess for
the phenology and population size. The initial system is 
is when both habitats have the same plant phenology. Our
initial guess of the mismatch in each trait is
```
> mm0=[-8.1471e+00     -8.0800e+00];
```
meaning the early habitat will hatch 8.1 days earlier than
the optimum, and the late habitat 8.08 days earlier than the
optium.

Our initial guess for the population size is
```
> n0=[1.8549e+01      7.4211e+01];
```
which means 18.5 individuals in the early habitat (remember
it starts with 20 territories, as defined in ```parV```) and 
74.2 individuals in the late habitat (it starts with 80
territories).

We can then use the function ```varypar.m``` to record the
singular strategy at each point in our parameter range
```
> varypar(p,'K',parV,1,mm0,n0)
```
The function in ```varypar.m``` produces an output file
called (in this case) ```outputK_1.dat```. 

## Reproduce the stability checks corresponding to Figure 2a

I have written a script ```check_stab_K_script.m``` that
will check the stability of each singular strategy found in the example
above for varying ```K``` in the early habitat.
