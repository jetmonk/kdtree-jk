# KDTREE-JK is a package for building efficent KD-Trees in Common Lisp

Data are stored in block arrays, for minimal (usually zero) consing
and allocation overheads.

Balancing is done in-place with Wirth's median partition method, so
should be fast and memory efficient.

By default, all floats are double precision, but this can be
changed in kdtree-jk-structs.lisp in the line
````
  (deftype kd-float ()
    'double-float)
````

## Main functions 

`(build-kdtree ndim :npoints 10)`

   Build a KDTREE of ndim dimensions, with initial allocation of
   npoints.  As points are added, the arrays are expanded by a factor of
   \*kdtree-expansion-factor\*, by default 1.5.  Thus it is nice but not
   crucial to  make npoints large enough to begin with.


`(insert-vector kdtree vec object &key defer)`

  Insert a vector V corresponding to OBJECT into NDIM KDTREE,
  returning the index where it ended up.

  If DEFER keyword is set, then the tree-insertion is deferred to save
  time, but it will be necessary to run BALANCE-KDTREE before
  use. This is useful for creating a balanced tree from unbalanced or
  nonrandom data.  Searching a tree that has not been balanced but
  that was inserted with DEFER generates an error.

`(insert-2d kdtree x y   object &key (vec nil) (defer nil))`

`(insert-3d kdtree x y z object &key (vec nil) (defer nil))`

  2d and 3d convenience functions.  VEC is a optional vector of the
  correct kd-float type to reduce consing, and DEFER is as above.

`(balance-kdtree kdtree)`

  Balance KDTREE.  Useful if data were inserted with :DEFER, or if
  data were non-random or unbalanced.  This is efficient, consing
  only one integer vector of the size of the tree.


`(kd-search-in-radius kdtree v radius &key (kdresult nil) (sort nil))`

  Search KDTREE for points within RADIUS of float vector V.  If
  keyword :KDRESULT is NIL, then a KDRESULT (see below) will be
  created.  Otherwise, given KDRESULT will be used to return the
  result.

  If SORT is true, then the KDRESULT is sorted by increasing distance
  from point given.

`(kd-search-in-box kdtree bbox &key (kdresult nil))`

  Search KDTREE for points inside bounding box BBOX (of type BBOX, see
  below).  If keyword :KDRESULT is NIL, then a KDRESULT will be
  created.  Otherwise, given KDRESULT will be used to return the
  result.

`(kd-find-k-nearest kdtree v k-nearest &key (rstart nil) (kdresult nil))`

  Find k-nearest objects in kdtree, by doing increasing radial
  searches. The returned KDRESULT will contain at least K-NEAREST
  points, sorted by radius.  It us up to the user to truncate the
  result to the first K-NEAREST.

  RSTART is an optional starting radius size.  By default, it it set
  using the volume fraction expected to be occupied by K-NEAREST points.

  Returns the number of search iterations as the second value.

`(kdtree-idepth kdtree)`

   The maximum depth of a KDTREE.  Should not be much more than
   (log npoints 2) if KDTREE is reasonably balanced.

`(kdtree-avg-depth kdtree)`

   Float representing the average depth of a node.a


## KDRESULT 

`(build-kdresult &key (n 20))`

  Create a KDRESULT with 20 slots; it will be expanded as necessary.

`(kdresult-n kdresult)`  -  The number of search results returned
`(kdresult-obj-vec kdresult)` - a vector of the OBJECTS returned
`(kdresult-index-vec kdresult)` - a vector of the indices in the kdtree
`(kdresult-dist-vec kdresult)` - a vector of distances from the search
                               vector in (kd-search-in-radius ...)
  
## BBOX (Bounding Box) 

`(build-bbox '(x1 y1 z1 ..) '(x2 y2 z2 ..))`

  Build a BBOX extending from x1..x2, y1..y2, etc, that is used
  to constrain a box search.  There is also a diagnostic BBOX inside
  a KDTREE to identify its covered space.


## Testing functions 

kdtree-jk-tests.lisp contains functions to compare box and radius
searches with 


