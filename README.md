# ComputationalPhysics2

 Git repository for the course "Computational Physics", master degree in Physics at UniTN.

## Course website

Click this [link](https://didatticaonline.unitn.it/dol/course/view.php?id=34003) to open the course web page on moodle.

## Repository contents

> The folder `projects` contains the various exercises, each in its own subfolder.
>
> The folder `myfunc` consists in a library with some functions to be used in all the exercises.
>
> The folder `reports` contains the reports of the exercises.

---

## Instructions to run the code

Open the terminal in the folder containing this README, and execute the following code:

```C
mkdir build
cd build
cmake ..
make -j
```

Then execute the instructions for the single projects, always remaining in the `build` folder.

---

## Project 0: test. Currently: Simple solving 1D SE with numerov

```C
make run-test
```

---

## Project 1: scattering H-Kr

To run the code for the scattering, type:

```C
make run-ex1-scatt
```

Instead to run the code for the bound states, type:

```C
make run-ex1-bound
```

To test the routine for calculating the Bessel functions type:

```C
make run-ex1-bessel_test
```

...
