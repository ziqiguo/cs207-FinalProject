# CS207 G8 Final Project: Chemical Kinetics

[![Build Status](https://travis-ci.org/G8-CS207F17/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/G8-CS207F17/cs207-FinalProject)
[![Coverage Status](https://coveralls.io/repos/github/G8-CS207F17/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/G8-CS207F17/cs207-FinalProject?branch=master)


Introduction
------------
Our package computes the progress rates of of a system chemical reactions. The package extracts parameters for each chemical reaction from an `xml` file. We also provide functions to calculate three types of reaction rate coefficients (Constant, Arrhenius and Modified Arrhenius) as well as progress rate for a reaction set. Based on the reaction rate coefficients and progress rates, the reaction rate of a system of reactions can be obtained.

Our package can also handle an arbitrary number of species and reactions. It is compatible with elementary and irreversible reactions, and can be easily extended to reversible or non-elementary reactions.




Installation
------------
Download package distribution from https://github.com/G8-CS207F17/cs207-FinalProject.

Packages needed: `numpy`, `xml`, `pytest` (`pip install PACKAGE_NAME`)

Run unit test suites: run `pytest --cov --cov-report term-missing` in terminal

run doctest suites: run `pytest --doctest-modules --cov --cov-report term-missing` in terminal



Basic Usage and Examples
------------------------
Provide a few examples on using your software in some common situations.  You may want to show how the code works with a small set of reactions.
