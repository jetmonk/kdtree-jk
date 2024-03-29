

(in-package kdtree-jk)


(defun build-bbox (min-bounds max-bounds)
  "Build a Bounding Box (BBOX).  MIN-VEC and MAX-VEC are sequences
representing the bounds."
  (declare (type sequence min-bounds max-bounds))
  (when (not (= (length min-bounds) (length max-bounds)))
    (error "MIN-BOUNDS and MAX-BOUNDS have different lengths"))
  (let* ((ndim (length min-bounds))
	 (vmin (make-kdfltvec ndim))
	 (vmax (make-kdfltvec ndim)))
    (loop for i below ndim
	  for xmin = (nth i min-bounds)
	  for xmax = (nth i max-bounds)
	  if (<= xmax xmin)
	    do (error "MIN-BOUNDS=~A is not less than MAX-BOUNDS=~A in each dimension."
		      min-bounds max-bounds)
	  else
	    do (setf (aref vmin i) (coerce  xmin 'kd-float))
	       (setf (aref vmax i) (coerce  xmax 'kd-float)))
      (make-bbox
       :ndim ndim
       :rmin-vec vmin
       :rmax-vec vmax)))

(defun build-kdresult (&key (n 20))
  (make-kdresult
   :n 0
   :obj-vec (make-objvec n)
   :index-vec (make-indexvec n)
   :dist-vec (make-kdfltvec n)))

(defun %expand-kdresult (kdresult)
  (declare (type kdresult kdresult)
	   (optimize speed))
  (let* ((n (kdresult-n kdresult))
	 (n-new (* 2 n))
	 (obj-vec-new (make-objvec n-new))
	 (index-vec-new (make-indexvec n-new))
	 (dist-vec-new (make-kdfltvec n-new)))
    (declare (type objvec obj-vec-new)
	     (type indexvec index-vec-new)
	     (type kdfltvec dist-vec-new))
    (loop for i of-type fixnum below n
	  do (setf (aref obj-vec-new i) (aref (kdresult-obj-vec kdresult) i))
	     (setf (aref index-vec-new i) (aref (kdresult-index-vec kdresult) i))
	     (setf (aref dist-vec-new i) (aref (kdresult-dist-vec kdresult) i)))
    (setf (kdresult-obj-vec kdresult) obj-vec-new)
    (setf (kdresult-index-vec kdresult) index-vec-new)
    (setf (kdresult-dist-vec kdresult) dist-vec-new)
    kdresult))
	     

;; insert index K from KDTREE into KDRESULT, at distance DIST
;; if OVERWRITE is set, then don't increment the counter so that
;; this point overwrites the last one
(declaim (inline %insert-point-in-kdresult))
(defun %insert-point-in-kdresult (kdtree kdresult k dist2 overwrite)
  (declare (type kdtree kdtree)
	   (type kdresult kdresult)
	   (type kd-float dist2)
	   (optimize (speed 3) (safety 1)))
  (when (= (kdresult-n kdresult)
	    (length (kdresult-obj-vec kdresult)))
    (%expand-kdresult kdresult))
  (let ((n (kdresult-n kdresult))
	(j (aref (kdtree-ir-vec kdtree) k))) ;; indirection
    (setf (aref (kdresult-obj-vec kdresult) n)
	  (aref (kdtree-obj-vec kdtree) j))
    (setf (aref (kdresult-index-vec kdresult) n) k)
    (setf (aref (kdresult-dist-vec kdresult) n) (sqrt (the (real 0) dist2)))
    (when (not overwrite)
      (incf (kdresult-n kdresult)))))

  

;; FIXME - would be nice to have specialized in-place heapsort without consing, but
;; GC time is negligible
(defun sort-kdresult-by-radius (kdresult)
  (when (> (kdresult-n kdresult) 1) ;; just 1 point might be a common result
    (let ((vec (loop with v = (make-array (kdresult-n kdresult))
		     for i below (kdresult-n kdresult)
		     do (setf (aref v i) (list (aref (kdresult-dist-vec kdresult) i)
					       (aref (kdresult-index-vec kdresult) i)
					       (aref (kdresult-obj-vec kdresult) i)))
		     finally (return v))))
      (setf vec (sort vec #'< :key #'first)) ;; sort by radius
      (loop for i below (kdresult-n kdresult)
	    for lst = (aref vec i)
	    do (setf (aref (kdresult-dist-vec kdresult)  i) (first lst)
		     (aref (kdresult-index-vec kdresult) i) (second lst)
		     (aref (kdresult-obj-vec kdresult)   i) (third lst))))))



(defun kd-search-in-radius (kdtree v radius &key (kdresult nil) (sort nil) (nearest-point nil) (action nil))
  "Search KDTREE for points within RADIUS of float vector V.  If
keyword :KDRESULT is NIL, then a KDRESULT will be created.  Otherwise,
given KDRESULT will be used to return the result.

If NEAREST-POINT is set, return just the nearest point.  This is
especially efficient because the internal search radius keeps
shrinking.

If SORT is true, then the KDRESULT is sorted by increasing distance
from point given.

If ACTION is supplied, than it is applied to each matching point, and
KDRESULT is ignored and set to NIL.  
It is called as (FUNCALL ACTION INODE) where INODE is the linear index of a node.

:ACTION (MAKE-DELETION-ACTION), for example, deletes any element returned by
search by setting its value to a special symbol 'DELETED-OBJECT."
  (declare (type kdtree   kdtree)
	   (type kdfltvec v)
	   (type (or null kdresult))
	   (type (or function symbol) action)
	   (type kd-float radius)
	   (optimize speed))
  (when (kdtree-needs-balancing kdtree)
    (error
     "KDTREE needs balancing with (BALANCE-KDTREE KDTREE) because insertions were performed with :DEFER"))
  (let* ((kdresult (if action nil (or kdresult (build-kdresult :n 10))))
	 (%action (if (and action (symbolp action)) (fdefinition action) action)) 
	 (ir-vec (kdtree-ir-vec kdtree))
	 (r-vec (kdtree-r-vec kdtree))
	 (obj-vec (kdtree-obj-vec kdtree))
	 (index-left-vec (kdtree-index-left-vec kdtree))
	 (index-right-vec (kdtree-index-right-vec kdtree))
	 (ndim (kdtree-ndim kdtree))
	 (nd-1 (1- ndim))
	 ;; search distance, which may shrink if we're finding just NEAREST-POINT
	 (dist2 (if (>= radius #.(sqrt +most-positive-float+))
		    +most-positive-float+ ;; avoid overflow for nearest-search case
		    (expt radius 2))))
    (declare (type kd-float dist2)
	     (fixnum ndim nd-1)
	     (type (or null function) %action)
	     (optimize (speed 3) (safety 1))) ;; safety=0 provides 14% speed boost in SBCL
    (when kdresult (setf (kdresult-n kdresult) 0))
    (labels ((%find-nearest (inode idim) ;; recursive function 
	       (when (= inode +end+) (return-from %find-nearest))
	       (let ((r2 (loop with j = (aref ir-vec inode) ;; indirection
			       for  idim below ndim
			       sum (expt (- (aref v idim) (aref r-vec j idim)) 2)
				 of-type kd-float))
		     (j (aref ir-vec inode))) ;; indirection
		 (declare (type kd-float r2))
		 
		 (when (<= r2 dist2)
		   (if action
		       (funcall (the function %action) kdtree inode)
		       (when (not (eq (aref obj-vec j) 'deleted-object))
			 (%insert-point-in-kdresult
			  kdtree kdresult inode r2
			  nearest-point) ;; overwrite last point if nearest-point=T
			 (when nearest-point
			   ;; make the search radius smaller than the best point
			   (setf dist2 r2)))))
		 
		 ;; search the side of the tree that our point is on
		 (let* ((dx (- (aref v idim) (aref r-vec j idim))))
		   (declare (type kd-float dx))
		   
		   (%find-nearest 
		    (if (<= dx 0)
			(aref index-left-vec inode)
			(aref index-right-vec inode))
		    (if (= idim nd-1) 0 (1+ idim)))

		   ;; if our point is within radius of the splitting x,y,z,
		   ;; then we have to search the other side of the tree too
		   (when (< (* dx dx) dist2)
		     (%find-nearest
		      (if (<= dx 0)
			  (aref index-right-vec inode)
			  (aref index-left-vec inode))
		      (if (= idim nd-1) 0 (1+ idim))))))))

      (%find-nearest 0 0) ;; start at first node
      (when (and sort (not nearest-point)) (sort-kdresult-by-radius kdresult))
      (when nearest-point (setf (kdresult-n kdresult) 1))
      kdresult)))

(defun kd-find-nearest-point (kdtree v &key (kdresult nil))
  "Find the single nearest point to vector V in KDTREE. 
A wrapper for:
  (KD-SEARCH-IN-RADIUS KDTREE V HUGE-RADIUS :KDRESULT KDRESULT NEAREST-POINT T)

The ACTION keyword is not supported because the result is not known
before the end of the search."
  (kd-search-in-radius kdtree v +most-positive-float+
		       :kdresult  (or kdresult (build-kdresult :n 1))
		       :nearest-point t))
  

(defun kd-find-k-nearest (kdtree v k-nearest &key (rstart nil) (kdresult nil))
  "Find k-nearest objects in KDTREE, by doing increasing radial
searches. The returned KDRESULT will contain at least K-NEAREST points, sorted
by radius.  It us up to the user to truncate the result to the first
K-NEAREST.

RSTART is an optional starting radius size.  By default, it it set using 
the volume fraction expected to be occupied by K-NEAREST points.

Returns the number of search iterations as the second value.

The ACTION keyword is not supported because the result is not know
before the end of the search."
  (declare (type kdtree   kdtree)
	   (type kdfltvec v)
	   (type (integer 0)  k-nearest)
	   (type (or null kdresult)))
  (when (> k-nearest (kdtree-npoints kdtree))
    (error "Cannot obtain ~D K-nearest points from a KDTREE with ~D elements."
	   k-nearest (kdtree-npoints kdtree)))
  (let* ((kdresult (or kdresult (build-kdresult :n (* 2 k-nearest))))
	 ;; fraction of points being retrieved
	 (frac (/ k-nearest (* 1.0 (kdtree-npoints kdtree))))
	 (d-largest
	   (loop with bbox = (kdtree-bbox kdtree)
		 for r1 across (bbox-rmin-vec bbox)
		 for r2 across (bbox-rmax-vec bbox)
		 maximizing (- r2 r1) into maxval
		 finally
		    (if (plusp maxval) 
			(return maxval)  ;; our other choice
			(error "The points in KDTREE cover no volume!"))))
	 (rsearch0
	   ;; set the initial search radius to be roughly the right
	   ;; volume fraction
	   (coerce (or rstart
		       (* d-largest 0.4 (expt frac 0.333)))
		   'kd-float)))

    (loop for rsearch = rsearch0 then (* 2 rsearch)
	  for itries from 1
	  for kdr = (kd-search-in-radius kdtree v rsearch
					      :kdresult kdresult :sort nil)
	  when (>= (kdresult-n kdresult) k-nearest)
	    do (sort-kdresult-by-radius kdr)
	       (return (values kdr itries)))))
	       
			   


(defun kd-search-in-bbox (kdtree bbox &key (kdresult nil) (action nil))
  (declare (type kdtree   kdtree)
	   (type (or null kdresult))
	   (type (or function symbol) action)
	   (type bbox bbox))
    "Search KDTREE for points inside bounding box BBOX (of type BBOX).  If keyword :KDRESULT
is NIL, then a KDRESULT will be created.  Otherwise, given KDRESULT will be used to
return the result.

If ACTION is supplied, than it is applied to each matching point, and
KDRESULT is ignored and set to NIL.  
It is called as (FUNCALL ACTION INODE) where INODE is the linear index of a node.

:ACTION (MAKE-DELETION-ACTION), for example, deletes an element by setting its value
to a special symbol 'DELETED-OBJECT."
  (when (kdtree-needs-balancing kdtree)
    (error
     "KDTREE needs balancing with (BALANCE-KDTREE KDTREE) because insertions were performed with :DEFER"))
  (let* ((kdresult (if action nil (or kdresult (build-kdresult :n 10))))
	 (%action (if (and action (symbolp action)) (fdefinition action) action)) 
	 (bbrmin (bbox-rmin-vec bbox))
	 (bbrmax (bbox-rmax-vec bbox))
	 (ir-vec (kdtree-ir-vec kdtree))
	 (r-vec (kdtree-r-vec kdtree))
	 (obj-vec (kdtree-obj-vec kdtree))
	 (index-left-vec (kdtree-index-left-vec kdtree))
	 (index-right-vec (kdtree-index-right-vec kdtree))
	 (ndim (kdtree-ndim kdtree))
	 (nd-1 (1- ndim)))
    (declare (optimize (speed 3) (safety 1))
	     (type (or null function) %action))
    (when kdresult (setf (kdresult-n kdresult) 0))
    (when (not (= ndim (bbox-ndim bbox)))
      (error "Bounding Box BBOX has wrong number of dimensions for KDTREE."))
    (labels ((%find-inbox (inode idim) ;; recursive function 
	       (when (= inode +end+) (return-from %find-inbox))
	       (let* ((j (aref ir-vec inode)) ;; indirection
		      ;; precompute some common array references
		      (rthis (aref r-vec j idim)) ;; value of the this dim
		      (boxmin (aref bbrmin idim))
		      (boxmax (aref bbrmax idim))
		      (in-box ;; is this node in the box?
			(loop for jdim below ndim
			      when (not (<= (aref bbrmin jdim)
					    (aref r-vec j jdim)
					    (aref bbrmax jdim)))
				do (return nil)
			      finally (return t))))
		 (declare (type kd-float rthis boxmin boxmax))
		 (when in-box
		   (if action
		       (funcall (the function %action) kdtree inode)
		       (when (not (eq (aref obj-vec j) 'deleted-object))
			 (%insert-point-in-kdresult kdtree kdresult inode
						    (coerce 0 'kd-float) nil))))

		 
		 ;; dx: which side of box is this node on?
		 (let* ((dx (- (/ (+ boxmin boxmax) 2)
			       rthis))
			;; search the other side on overlap of the box with the node
			(search-other-box
			  (or
			   (and (<= dx 0) (> boxmax rthis))
			   (and (>= dx 0) (< boxmin rthis)))))
		   (declare (type kd-float dx))

		   (%find-inbox
		    (if (<= dx 0)
			(aref index-left-vec inode)
			(aref index-right-vec inode))
		    (if (= idim nd-1) 0 (1+ idim)))

		   ;; if our point is within radius of the splitting x,y,z,
		   ;; then we have to search the other side of the tree too
		   (when search-other-box
		     (%find-inbox
		      (if (<= dx 0)
			  (aref index-right-vec inode)
			  (aref index-left-vec inode))
		      (if (= idim nd-1) 0 (1+ idim))))))))
	     
      (%find-inbox 0 0) ;; start at first node
      kdresult)))





