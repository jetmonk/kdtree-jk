#|
   run KD tree creation and search on random data and and random bboxes
|#


(defpackage kdtree-jk/tests
  (:use #:cl #:kdtree-jk)
  (:export
   #:run-radius-test-set
   #:run-bbox-test-set
   #:run-big-test-series
   #:wirth-partition-test
   ))

(in-package kdtree-jk/tests)

(defun make-data (&key (ndim 3) (npts 1000) (range-min -100d0) (range-max 100d0))
  (loop for i below npts
	collect (loop with v = (kdtree-jk::make-kdfltvec ndim)
		      for i below ndim
		      for x = (coerce (+ range-min (random (- range-max range-min)))
				      'kdtree-jk::kd-float)
		      do (setf (aref v i) x)
		      finally (return v))))

(defun build-test-kdtree (&key (ndim 3) (npts 1000) (range-min -100d0) (range-max 100d0)
			       (balance nil))
  (let ((data (make-data :ndim ndim :npts npts :range-min range-min :range-max range-max))
	(kdtree (build-kdtree ndim :npoints 10))) ;; make it expand itself
    (loop for vec in data
	  do (insert-vector kdtree vec vec :defer balance))
    (when balance (kdtree-jk::balance-kdtree kdtree))
    (values kdtree data)))

(defun %vec-dist (vec1 vec2)
  (loop for x1 across vec1 and x2 across vec2
	sum (expt (- x1 x2) 2)))

;; generate a hash of those vectors in vec-list that are within r of vec
(defun %brute-force-match-in-radius (vec vec-list r)
  (loop with hash = (make-hash-table :test 'eq)
	for vec2 in vec-list
	for radius2 = (%vec-dist vec vec2)
	when (<= radius2 (* r r))
	  do (setf (gethash vec2 hash) t)
	finally (return hash)))

(defun run-one-radius-test (&key (ndim 3) (npts 1000) (range-min -100d0) (range-max 100d0)
				 (radius 20) (balance nil))
  (multiple-value-bind (kdtree kd-vec-list)
      (build-test-kdtree :ndim ndim :npts npts :range-min range-min :range-max range-max :balance balance)
    (let ((other-data ;; the set we compare against
	    (make-data :ndim ndim :npts npts :range-min range-min :range-max range-max))
	  (total-matches 0))
      (loop for vec in other-data
	    for kdresult = (kd-search-in-radius kdtree vec (coerce radius 'kdtree-jk::kd-float))
	    for brute-hash = (%brute-force-match-in-radius vec kd-vec-list (* 1d0 radius))
	    do
	       (incf total-matches (kdresult-n kdresult))
	       (when (not (= (hash-table-count brute-hash)
			     (kdresult-n kdresult)))
		 (error "Number of matches returned by brute force is ~A but ~A by KD match"
			(hash-table-count brute-hash) (kdresult-n kdresult)))
	       (loop for i below (kdresult-n kdresult)
		     for obj across (kdresult-obj-vec kdresult)
		     when (not (gethash obj brute-hash))
		       do (error "Object ~A is not found in brute force match."
				 obj)))
      total-matches)))

(defun run-radius-test-set (&key (ndim 3) (npts 1000) (balance nil))
  (format t "If this does not produce an error, then a KD radial search produces the same
result as a brute force match.~%~%")
  (loop with ntot = 0
	for radius from 1.0 to 100.0
	do (incf ntot (run-one-radius-test :ndim ndim :npts npts :balance balance))
	finally
	   (format t "Successfully matched ~D total objects in agreement using KD and brute force~%"
		   ntot)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; is a vector in a bbox?
(defun %vec-in-bbox? (vec bbox)
  (loop for i below (bbox-ndim bbox)
	when (not (<= (aref (bbox-rmin-vec bbox) i)
		      (aref vec i)
		      (aref (bbox-rmax-vec bbox) i)))
	  do (return nil)
	finally (return t)))
		      

(defun %brute-force-match-in-bbox (bbox vec-list)
  (loop with hash = (make-hash-table :test 'eq)
	for vec2 in vec-list
	when (%vec-in-bbox? vec2 bbox)
	  do (setf (gethash vec2 hash) t)
	finally (return hash)))

(defun %build-random-bbox (ndim range-min range-max bbox-scale)
  (let* ((size (- range-max range-min))
	 (r (+ 0.1 (random 1.0))) ;; 0.1 because we don't want super-thin boxes
	 ;; start the box in some random region of [range-min,range-max]
	 (bbox-min-list (loop for i below ndim
			      collect (+ range-min (* size (random 0.8)))))
	 ;; end the box in the start, plus some random factor of bbox-scale
	 (bbox-max-list (loop for bmin in bbox-min-list
			      collect (+ bmin (* size r bbox-scale)))))
    (build-bbox bbox-min-list bbox-max-list)))
	 


;; this uses larger npoints because the boxes are generally smaller than large radii
(defun run-one-bbox-test (&key (ndim 3) (npts 100000) (range-min -100d0) (range-max 100d0)
			    (bbox-scale 0.01) (nbbox 1000) ;; size of region
			    (balance nil))
  (multiple-value-bind (kdtree kd-vec-list)
      (build-test-kdtree :ndim ndim :npts npts :range-min range-min :range-max range-max :balance balance)
    (let ((total-matches 0))
      (loop for j below nbbox 
	    for bbox = (%build-random-bbox ndim range-min range-max bbox-scale)
	    for kdresult = (kd-search-in-bbox  kdtree  bbox)
	    for brute-hash = (%brute-force-match-in-bbox bbox kd-vec-list)
	    do
	       (incf total-matches (kdresult-n kdresult))
	       (when (not (= (hash-table-count brute-hash)
			     (kdresult-n kdresult)))
		 (error "Number of matches returned by brute force is ~A but ~A by KD match"
			(hash-table-count brute-hash) (kdresult-n kdresult)))
	       (loop for i below (kdresult-n kdresult)
		     for obj across (kdresult-obj-vec kdresult)
		     when (not (gethash obj brute-hash))
		       do (error "Object ~A is not found in brute force match."
				 obj)))
      total-matches)))

  
(defun run-bbox-test-set (&key (ndim 3) (npts 10000) (nbbox 1000) ;; will take a long time to run
			       (balance nil))
  (format t "If this does not produce an error, then a KD BBOX search produces the same
result as a brute force match.~%~%")
  (loop
    with ntot = 0
    for bbox-scale from 0.5 to 1.0 by 0.05
    do (incf ntot (run-one-bbox-test :ndim ndim :npts npts :nbbox nbbox :balance balance))
    finally (format t "Successfully matched ~D total objects in agreement using KD and brute force~%"
		    ntot)))
  

;; this might be an overnight test
(defun run-big-test-series (&key (balance nil))
  (loop for ndim in '(2 3 4) ;; higher will produce no hits
	do (run-radius-test-set :ndim ndim :npts 10000 :balance balance))
  (loop for ndim in '(2 3 4) ;; higher will produce no hits
	do (run-bbox-test-set :ndim ndim :npts 100000 :nbbox 1000 :balance balance)))
  
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defun %verify-wirth-partition (irvec rvec idim j j1 j2)
  (declare
   (type kdtree-jk::indexvec irvec)
   (type kdtree-jk::kdfltarr rvec)
   (type fixnum idim j j1 j2)
   (optimize speed))
  (loop with nbelow of-type fixnum = 0
	with nabove of-type fixnum = 0
	with rj = (aref rvec (aref irvec j) idim)
	for i from j1 to j2
	for ri  = (aref rvec (aref irvec i) idim)
	do (cond ((< i j)
		  (if (> ri rj)
		      (error "Not ordered correctly (below)"))
		  (incf nbelow))
		 ((> i j)
		  (if (< ri rj)
		      (error "Not ordered correctly (above)"))
		  (incf nabove)))
	finally
	   (if (> (abs (- nabove nbelow)) 1)
	       (error "Not a median because |NABOVE-NBELOW|>1"))))
  


(defun wirth-partition-test (&key (niter 1000) (ndim 3) (n 1000))
  "Test to see if partitioning scheme to find median and subdivide inputs gives errors"
  (loop for i below niter
	for j1 = (random (ash n -1))
	for j2 = (min (1- n) (+ j1 (random n)))
	for kd = (build-kdtree ndim :npoints n)
	for idim = (random ndim)
	do 
	   (loop with v = (kdtree-jk::make-kdfltvec ndim)
		 for i below 1000
		 do 
		    (loop for kdim below ndim do
		      (setf (aref v kdim) (coerce (random 10d0) 'kdtree-jk::kd-float)))
		    (kdtree-jk::insert-vector kd v t))
	   (let ((j (kdtree-jk::partition-ir-vec-around-median kd j1 j2 idim)))
	     (%verify-wirth-partition (kdtree-jk::kdtree-ir-vec kd) (kdtree-jk::kdtree-r-vec kd)
				      idim
				      j j1 j2))))
	
