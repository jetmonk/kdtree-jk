# KDTREE-JK is a package for building efficent KD-Trees in Common Lisp

A KD-Tree is a data structure for searching N-dimensional space in
log(N) time by partitioning a data set iteratively across dimensions.

In KDTREE-JK, data are stored in block arrays instead of linked
structures, for minimal (usually zero) consing and allocation
overheads.  This is in contrast to a canonical KD-Tree implementation
that uses linked lists, which wastes memory on boxing.

Balancing is done in-place with Wirth's median partition method, so it
should be fast and memory efficient.  Insertion of one object and
searching is O(log(N)), and re-partitioning is O(N[log(N)]^2)

By default, all floats are double precision, but this can be
changed in kdtree-jk-structs.lisp in the line
````
  (deftype kd-float ()
    'double-float)
````

## Documentation

See `doc/reference.txt`

See doc/reference-latlon.txt for documentation to a second package,
    KDTREE-JK/LATLON, to search in longitude,latitude coordinates
    by converting to 3d points on a unit sphere.

kdtree-jk is written in pure Common Lisp, with no dependencies.

Submit bug reports at https://github.com/jetmonk/kdtree-jk/issues