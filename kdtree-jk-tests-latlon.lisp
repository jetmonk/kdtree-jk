
(in-package kdtree-jk/tests)


(defun xyz->lon-lat (x y z)
  (let* ((r (sqrt (+ (expt x 2) (expt y 2) (expt z 2))))
	 (r% (if (zerop r) 1d0 r))
	 (lat (asin (/ z r%)))
	 (lon (atan y x))
	 (180/pi (/ 180 pi)))
    (declare (type double-float r r% lat lon)
               (dynamic-extent r% lat lon))
    (setf lat (* lat 180/pi)
	  lon (* lon 180/pi))
      (values (coerce lon 'kdtree-jk:kd-float)
	      (coerce lat 'kdtree-jk:kd-float))))

(defun generate-random-lon-lat ()
  (xyz->lon-lat
   (+ -1d0 (* 2d0 (random 1d0)))
   (+ -1d0 (* 2d0 (random 1d0)))
   (+ -1d0 (* 2d0 (random 1d0))))) 

;; uses formula for angle more accurate than acos(dotproduct)
(defun sp-angle (lon1 lat1 lon2 lat2)
  (let ((v1 (kdtree-jk::make-kdfltvec 3))
	(v2 (kdtree-jk::make-kdfltvec 3)))
    (kdtree-jk/latlon::lon-lat-to-xyz-vec lon1 lat1 v1)
    (kdtree-jk/latlon::lon-lat-to-xyz-vec lon2 lat2 v2)
    (let* ((x1 (aref v1 0))
	   (y1 (aref v1 1))
	   (z1 (aref v1 2))
	   (x2 (aref v2 0))
	   (y2 (aref v2 1))
	   (z2 (aref v2 2))
	   (delta/2 (min 1d0 ;; avoid (asin 1.00000000001)
			 (* 0.5d0 (sqrt (+ (expt (- x1 x2) 2)
					   (expt (- y1 y2) 2)
					   (expt (- z1 z2) 2)))))))
      (* 2d0 (/ 180 pi) (asin (the (double-float 0d0 1d0) delta/2))))))
      


(defun make-sim-data-latlon (&key  (npts 1000))
  (loop for i below npts
	collect (multiple-value-bind (lon lat)
		    (generate-random-lon-lat)
		  (vector lon lat))))



(defun %brute-force-match-in-latlon (vec vec-list angle)
  (loop with hash = (make-hash-table :test 'eq)
	for vec2 in vec-list
	for dist = (sp-angle (aref vec 0) (aref vec 1)
			     (aref vec2 0) (aref vec2 1))
	when (<= dist angle)
	  do (setf (gethash vec2 hash) (sqrt dist))
	finally (return hash)))


(defun build-test-latlon-kdtree (&key  (npts 1000)
				   (balance nil))
  (let ((data (make-sim-data-latlon  :npts npts))
	(kdtree (build-kdtree-latlon :npoints 10))) ;; make it expand itself
    (loop for vec in data
	  do (insert-latlon  kdtree  (aref vec 0) (aref vec 1) vec :defer balance))
    (when balance (kdtree-jk::balance-kdtree kdtree))
    (values kdtree data)))

(defun run-one-latlon-test (&key  (npts 1000)
				 (angle 10.0) (balance nil))
  (multiple-value-bind (kdtree kd-vec-list)
      (build-test-latlon-kdtree :npts npts  :balance balance)
    (let ((other-data ;; the set we compare against
	    (make-sim-data-latlon  :npts npts))
	  (total-matches 0))
      (loop for vec in other-data
	    for iter from 0
	    for kdresult = (kd-latlon-search-in-angle  kdtree 
						       (aref vec 0) (aref vec 1)
						       (coerce angle 'kdtree-jk::kd-float))
	    for brute-hash = (%brute-force-match-in-latlon vec kd-vec-list (* 1d0 angle))
	    do
	       (incf total-matches (kdresult-n kdresult))
	       (when (not (= (hash-table-count brute-hash)
			     (kdresult-n kdresult)))
		 (error "Number of matches returned by brute force is ~A but ~A by KD match in iter ~A"
			(hash-table-count brute-hash) (kdresult-n kdresult) iter))
	       (loop for i below (kdresult-n kdresult)
		     for obj across (kdresult-obj-vec kdresult)
		     for dist across (kdresult-dist-vec kdresult)
		     when (not (gethash obj brute-hash))
		       do (error "Object ~A at dist=~A is not found in brute force match."
				 obj dist)))
      total-matches)))

(defun run-latlon-test-set (&key (npts 1000) (balance nil))
  (format t "If this does not produce an error, then a KD lat-lon
search produces the same result as a brute force match.  It is
possible that trig roundoff might very rarely produce an apparent
failure.~%~%")
  (loop with ntot = 0
	for dist from 1.0 to 90.0 by 2.0
	do (incf ntot (run-one-latlon-test :npts npts :balance balance))
	finally
	   (format t "Successfully matched ~D total objects in agreement using KD and brute force~%"
		   ntot)))



