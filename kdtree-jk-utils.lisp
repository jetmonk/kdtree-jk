
(in-package kdtree-jk)

#|

Reference:
   Author: Wirth, Niklaus
   Title: Algorithms + data structures = programs
   Publisher: Englewood Cliffs: Prentice-Hall, 1976
   Physical description: 366 p.
   Series: Prentice-Hall Series in Automatic Computation

|#



#|

Given a KDTREE, this will rearrange the elements of elements of
KDTREE-IR-VEC between inclusive indices [J1,J2] so that J is the median
of the IDIM dimension of RVEC=KDTREE-R-VEC, and so for all J=J1..J2 the following
holds:

RVEC[IRVEC[J1:J-1],IDIM]  <=  RVEC[IRVEC[J],IDIM] <= RVEC[IRVEC[J+1:J2],IDIM]

Elements outside the range [J1:J2] are not touched.

Returns J and rearrnges KD-IR-VEC between J1 and J2 inclusive.

The purpose of this function is to partition the data in a way suitable for
balancing the tree recursively, so that J1:J-1 go on left side, and J+1:J2
go on right side, and J goes in the parent node.
|#

(defun partition-ir-vec-around-median (kdtree j1 j2 idim)

  (let ((irvec (kdtree-ir-vec kdtree))
	(rvec (kdtree-r-vec kdtree))) ;; 2d array []

    (labels ((wirth-index-of-kth-smallest (k nv1 nv2)
	       ;; kth smallest number with indices in [nv1,nv3], but also
	       ;; if the final index is kk, then this also
	       ;; ensures that v[j<kk]<=v[kk]<=v[j>kk]
	       ;; 
	       (declare (type fixnum k nv1 nv2)
			(optimize speed))
	       (let ((i 0) (j 0) (l nv1) (m 0) (x (coerce 0 'kd-float))
		     (kk (+ k nv1))
		     (n (1+ nv2)))
		 (declare (type fixnum i j l m n kk)
			  (type kd-float x))
		 (when (>= k n) 
		   (error "k=~D is too big; must be <~D" k n))
		 (setf m (1- n))
		 (loop
		   while (< l m)
		   do
		      (setf x (aref rvec (aref irvec kk) idim))
		      (setf i l)
		      (setf j m)
		      (loop
			do
			   (loop while (< (aref rvec (aref irvec i) idim) x) do (incf i))
			   (loop while (< x (aref rvec (aref irvec j) idim)) do (decf j))
			   (when (<= i j)
			     (rotatef (aref irvec i) (aref irvec j))
			     (incf i)
			     (decf j))
			while 
			(<= i j))
		      (when (< j kk) (setf l i))
		      (when (< kk i) (setf m j)))
		 ;;
		 kk)))
      (let ((j (wirth-index-of-kth-smallest  (ash (- j2 j1) -1) j1 j2)))
	;; (%verify-wirth-partition irvec rvec idim j j1 j2) ;; see kdtree-jk-tests.lisp
	j))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; remove the 'deleted-objects from kd-tree, compacting the
;; kdtree-r-vec and kdtree-obj-vec arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun %remove-deleted-objects (kdtree)
  (declare (type kdtree kdtree)
	   (optimize speed))
  (setf (kdtree-needs-balancing kdtree) t) ;; this function messes up the tree
  (let ((np (kdtree-npoints kdtree))
	(npnew 0)
	(rvec (kdtree-r-vec kdtree))
	(fix-bbox nil) ;; need to fix bounding box to reflect removed points
	(objvec (kdtree-obj-vec kdtree)))
    (declare (type index npnew))
    (loop with jtarg of-type index = 0
	  for  jsrc of-type index below np
	  if (not (eq (aref objvec jsrc) 'deleted-object))
	    do (incf npnew)
	       (when (not (= jsrc jtarg))
		 (setf (aref objvec jtarg)
		       (aref objvec jsrc))
		 (loop for idim below (kdtree-ndim kdtree)
		       do (setf (aref rvec jtarg idim)
				(aref rvec jsrc idim))))
	       (incf jtarg)
	  else
	    do (setf fix-bbox t))
    (setf (kdtree-npoints kdtree) npnew)

    ;; not necessary, but it will make kdtree easier to examine visually
    ;; and shorter to dump
    (loop with zero = (coerce 0 'kd-float)
	  for i from npnew below (kdtree-npoints-total kdtree)
	  do (setf (aref objvec i) nil)
	     (setf (aref (kdtree-ir-vec kdtree) i) +end+)
	     (setf (aref (kdtree-index-left-vec kdtree) i) +end+)
	     (setf (aref (kdtree-index-right-vec kdtree) i) +end+)
	     (setf (aref (kdtree-flag-vec kdtree) i) +empty+)
	     (setf (aref (kdtree-sd-vec kdtree) i) 0)
	     (loop for k below (kdtree-ndim kdtree)
		   do (setf (aref rvec i k) zero)))

    (when fix-bbox
      (loop with rmin-vec = (bbox-rmin-vec (kdtree-bbox kdtree))
	    with rmax-vec = (bbox-rmax-vec (kdtree-bbox kdtree))
	      initially (loop for i below (kdtree-ndim kdtree)
			      do (setf (aref rmin-vec i) (aref rvec 0 i))
				 (setf (aref rmax-vec i) (aref rvec 0 i)))
	    for j from 1 below npnew
	    do
	       (loop
		 for i below (kdtree-ndim kdtree)
		 for rji of-type kd-float = (aref rvec j i)
		 do (setf (aref rmin-vec i) (min rji (aref rmin-vec i)))
		    (setf (aref rmax-vec i) (max rji (aref rmax-vec i))))))
    t))
	  
	
