# Computational-physics library

Computational physics routines in C.

This repository contains some C routines for computational physics such as numerical integration, methods for differential equations and linear algebra wrappers.

The CMake structure creates a **static library** `myCP-lib` that you can install to your computer following the instructions in this file.

---

## Install

To install the code of this library locally, first download this repository from github.

After downloading, open a terminal and move to your local copy of the repository.

From the main directory, execute

```C
mkdir build
cd build
cmake ..
make -j
```

This will compile the library and put the executables in the directory `build`.

After compiling, execute the following command to install the library:

```C
make install
```

After doing that, the library will be installed in the location

```C
$HOME/.local/
```

After this procedure, you can include with brackets the header files of the `myCP-lib`.

---

## Linking (to verify)

Make sure to link the library when you compile, if you are using the `myCP-lib`. You can do this by compiling with the following flags:

```C
gcc <c_file> -o <executable_name> -I$HOME/.local/include -L$HOME/.local/lib -libmyCP-lib
```

---

## Content

The source code (i.e. the .c files) is in the directory `src`.

The header files are in the directory `include`.

---

## Documentation (also this is to verify)

The source code is documented with Doxygen. To generate the documentation, make sure to have Doxygen installed and then use `cmake` with the documentation option:

```C
cmake .. -DDOCS=ON
make -j
```

After that, you can find documentation files in the directory `build/docs`.
