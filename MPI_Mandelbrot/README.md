- Grid boundaries, the value which determines the end of point testing, and maximum number of iterations per point are defined in module _GlobalVars.f90_
- Throughout the program each point in a grid is referenced mainly through its global index which is recalculated into (i,j) poistion on the grid with function `index2point(...)`

- The logic of the program is the following:
  - Before parallelizing, root process (#0 in the program) scans the grid with a coarse step and fills array _grid_c_ for # of iterations per each point
  - After this, root process goes through a fine grid and calculates total weight for the grid. Total weight is the sum of weights for all points. And the wight of a given point equals to the number of iterations for the nearest point on the coarse grid (taken from array _grid_c_)
  - Total weight is divided by the number of processes to obtain "weight per process". Root process goes through the fine grid one more time and determines the limits of domains that must be attributed to each process by comparing the current weight of n points with the "weight per process". As soon as it exceeds, the current point determines the boundary for a process.
  - Finally, each process scans its domain and determines the number of iterations per each point. All data from processes is merged to form a single array (via MPI_GATHERV) which is passed to the function performing output to file
  - Logging (time of execution, number of iterations per each process)

- This division of tasks among process is fairly uniform: e.g. for 16 processes, step of coarse grid = 0.01 and step of fine grid = 0.001, the difference in the number of operations performed by each process is ~1 % (that can be seen from the _program_output.dat_)
- To see, how the required time changes, script _script.sh_ runs mpirun several times with different # of processes

- The [script](https://colab.research.google.com/drive/1wOJS__zrkOwJ6hAYOHo74vvqNVbHQov7?usp=sharing) in python to draw Mandelbrot set
