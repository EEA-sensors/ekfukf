EKF/UKF Toolbox for Matlab and GNU Octave
==

[Simo Särkkä](http://users.aalto.fi/~ssarkka/), Jouni Hartikainen, and [Arno Solin](http://arno.solin.fi)


Introduction
--
EKF/UKF is an optimal filtering toolbox for Matlab. Optimal filtering is a frequently used term for a process, in which the state of a dynamic system is estimated through noisy and indirect measurements. This toolbox mainly consists of Kalman filters and smoothers, which are the most common methods used in stochastic state-space estimation. The purpose of the toolbox is not to provide a highly optimized software package, but instead to provide a simple framework for building proof-of-concept implementations of optimal filters and smoothers to be used in practical applications.

Most of the code has been written by Prof. Simo Särkkä. Later Dr. Jouni Hartikainen and Arno Solin documented and extended it with new filters and smoothers as well as simulated examples.


Download and Installation Guide
--

Matlab
---

The software consists of Matlab m-files.
Clone or download the latest version and make sure the toolbox directory is included in
your Matlab path by `addpath` *path to ekfukf*.

GNU Octave
---

1. Run the following command

        make dist

    This will create a  file `ekfukf-<version>.tar.gz`

2. Install the package.
    Start octave and run

        pkg install ekfukf-<version>.tar.gz


Documentation
--
The documentation demonstrates the use of software as well as state-space estimation with Kalman filters in general. The purpose is not to give a complete guide to the subject, but to discuss the implementation and properties of Kalman filters.

* Documentation available on GitHub

The methods that are discussed in the current documentation are:

* Kalman filters and smoothers
* Extended Kalman filters and smoothers
* Unscented Kalman filters and smoothers
* Gauss-Hermite Kalman filters and smoothers
* Cubature Kalman filters and smoothers
* Interacting Multiple Model (IMM) filters and smoothers

Useful background information on the methods can also be found in the book:

* Simo Särkkä (2013). *Bayesian Filtering and Smoothing*. Cambridge University Press. Physical books available from [Cambridge University Press](http://www.cambridge.org/sarkka) and [a free electronic version online](http://users.aalto.fi/~ssarkka/).


Demos
--
There are a number of demonstration programs for the provided filters and smoothers. The code and a short introduction to them are given below. All of the demonstration programs are discussed in the documentation.

Demonstration programs for linear state-space models:

* 2D CWPA-model, `kf_cpwa_demo`

Demonstration programs for non-linear state-space models:

* Tracking a random sine signal, `ekf_sine_demo`
* UNGM-model, `ungm_demo`
* Bearings only tracking, `bot_demo`
* Reentry vehicle tracking, `reentry_demo`

Demonstration programs for multiple model systems:

* Tracking a target with simple manouvers, `imm_demo`
* Coordinated turn model, `ct_demo`
* Bearings only tracking of a manouvering target, `botm_demo`


License
--
This software is distributed under the GNU General Public License (version 2 or later); please refer to the file `LICENSE.txt`, included with the software, for details.
