# KDTREE-JK is a package for building efficent KD-Trees in Common Lisp

A KD-Tree is a data structure for searching N-dimensional space in
log(N) time by partiioning a data set iteratively across dimensions.

In KDTREE-JK, data are stored in block arrays instead of linked
structures, for minimal (usually zero) consing and allocation
overheads.

Balancing is done in-place with Wirth's median partition method, so it
should be fast and memory efficient.

By default, all floats are double precision, but this can be
changed in kdtree-jk-structs.lisp in the line
````
  (deftype kd-float ()
    'double-float)
````

## Documentation

See `doc/reference.txt`
