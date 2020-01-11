
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
	       (declare ;(type indexvecd  v)
			(type fixnum k nv1 nv2)
			(optimize speed))
	       (let ((i 0) (j 0) (l nv1) (m 0) (x 0.0)
		     (kk (+ k nv1))
		     (n (1+ nv2)))
		 (declare (type fixnum i j l m n kk)
			  (type single-float x))
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


