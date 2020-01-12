
(in-package kdtree-jk)

;; how much a kdtree expands in size when it runs out of memory
(defparameter *kdtree-expansion-factor* 1.5)

;; increase size of arrays in kdtree by FACTOR and copy data
;; if new-size is not NIL, then just use NEW-SIZE (can be
;; used to shrink kdtree)
(defun %expand-kdtree (kdtree &key
				(factor *kdtree-expansion-factor*)
				(new-size nil))
  (declare (type  (single-float 1.1 100.0) factor)
	   (type  (or null (unsigned-byte 32)) new-size)
           (optimize speed))
  (let* ((np (kdtree-npoints kdtree))
	 (npnew (or new-size (round (1+ (* factor np)))))
	 (ndim (kdtree-ndim kdtree))
	 (r-vec (kdtree-r-vec kdtree))
	 (flag-vec (kdtree-flag-vec kdtree))
	 (obj-vec (kdtree-obj-vec kdtree))
	 (index-left-vec (kdtree-index-left-vec kdtree))
	 (index-right-vec (kdtree-index-right-vec kdtree))
	 (ir-vec (kdtree-ir-vec kdtree))
	 (sd-vec (kdtree-sd-vec kdtree))
	 (r-vec-new (make-kdfltarr npnew ndim))
	 (flag-vec-new (make-flagvec npnew +empty+))
	 (obj-vec-new (make-objvec npnew nil))
	 (index-left-vec-new (make-indexvec npnew +end+))
	 (index-right-vec-new (make-indexvec npnew +end+))
	 (sd-vec-new (make-dimvec npnew))
	 (ir-vec-new (make-indexvec npnew +end+)))			      
    (declare (type indexvec index-left-vec index-left-vec-new
		   index-right-vec index-right-vec-new
		   ir-vec ir-vec-new)
	     (type kdfltarr r-vec r-vec-new)
	     (type objvec obj-vec obj-vec-new)
	     (type dimvec sd-vec sd-vec-new)
	     (type flagvec flag-vec flag-vec-new))
    (loop for i below np
	  do (setf (aref index-left-vec-new i)
		   (aref index-left-vec i))
	     (setf (aref index-right-vec-new i)
		   (aref index-right-vec i))
	     (setf (aref ir-vec-new i) (aref ir-vec i))
	     (setf (aref sd-vec-new i) (aref sd-vec i))
	     (setf (aref flag-vec-new i) (aref flag-vec i))
	     (setf (aref obj-vec-new i) (aref obj-vec i))
	     (loop for j below ndim
		   do (setf (aref r-vec-new i j) (aref r-vec i j))))
    (setf (kdtree-index-left-vec kdtree) index-left-vec-new)
    (setf (kdtree-index-right-vec kdtree) index-right-vec-new)
    (setf (kdtree-ir-vec kdtree) ir-vec-new)
    (setf (kdtree-sd-vec kdtree) sd-vec-new)
    (setf (kdtree-r-vec kdtree) r-vec-new)
    (setf (kdtree-flag-vec kdtree) flag-vec-new)
    (setf (kdtree-obj-vec kdtree)  obj-vec-new)
    (setf (kdtree-npoints-total kdtree) npnew)))


(defun kdtree-minimize-size (kdtree)
  "Make the arrays in KDTREE just large enough to hold the data."
  (%expand-kdtree kdtree :new-size (kdtree-npoints kdtree)))


;; jr is row index in r-vec, and k is index in
;; all other vectors; generally, they are the same, but
;; perhaps in a reshuffling they will be different
(declaim (inline %insert-into-kd)) 
(defun %insert-into-kd (v object kdtree jr k idepth idim)
  (declare (type kdfltvec v)
	   (type t object)
	   (type kdtree kdtree)
	   (type fixnum jr k idepth)
	   (type dimnum idepth)
	   (optimize speed))
  (incf (kdtree-npoints kdtree) 1)
  (setf (aref (kdtree-flag-vec kdtree) k) +filled+)
  (setf (aref (kdtree-sd-vec kdtree) k) idim)
  (setf (aref (kdtree-ir-vec kdtree) k) jr) ;; ir-vec points to r-vec
  (setf (aref (kdtree-obj-vec kdtree) k) object)
  (setf (kdtree-idepth kdtree) (max idepth (kdtree-idepth kdtree)))
  ;; compute the avg depth - not too accurate because of rounding
  (let ((f (/ 1d0 (kdtree-npoints kdtree)))) ;;fraction of points this point is
    (setf (kdtree-avg-depth kdtree)
	  (+ (* (- 1d0 f) (kdtree-avg-depth kdtree))
	     (* f idepth))))
  (let ((r-vec (kdtree-r-vec kdtree)))
    (loop with first-point = (= (kdtree-npoints kdtree) 1)
	  with rmin-vec = (bbox-rmin-vec (kdtree-bbox kdtree))
	  with rmax-vec = (bbox-rmax-vec (kdtree-bbox kdtree))	  
	  for i below (length v)
	  for r of-type kd-float = (aref v i)
	  do (setf (aref r-vec jr i) r)
	     (if first-point
		 (progn
		   (setf (aref rmin-vec i) r)
		   (setf (aref rmax-vec i) r))
		 (progn
		   (setf (aref rmin-vec i) (min r (aref rmin-vec i)))
		   (setf (aref rmax-vec i) (max r (aref rmax-vec i))))))))

;(defparameter *nops* 0)
;(declaim (type fixnum *nops*))

(defun insert-vector (kdtree v object &key (defer nil))
  "Insert a vector V corresponding to OBJECT into NDIM KDTREE,
returning the index where it ended up.

If DEFER keyword is set, then the tree-insertion is deferred to save
time, but it will be necessary to run BALANCE-KDTREE before use. This
is useful for creating a balanced tree from unbalanced or nonrandom
data."
  (declare (type kdfltvec v)
	   (type kdtree kdtree)
	   (optimize (speed 3) (safety 1)))
  (when (not (= (length v) (kdtree-ndim kdtree)))
    (error "Wrong number of dimensions"))
  ;; expand kdtree if necessary
  (when (= (kdtree-npoints-total kdtree)
	   (kdtree-npoints kdtree))
    (%expand-kdtree kdtree))

  ;; Is this a deferred addition? then just put it into the data
  ;; arrays, but don't walk the tree for the insert.  It will be
  ;; necessary to run BALANCE-KDTREE before accessing it.
  (when defer
    (setf (kdtree-needs-balancing kdtree) t)
    (%insert-into-kd v object kdtree
		     (kdtree-npoints kdtree)
		     (kdtree-npoints kdtree)
		     0 0)
    (return-from insert-vector nil))
  
  (loop with ndim = (kdtree-ndim kdtree)
	with k = 0 ;; current index in the tree
	with knf = (kdtree-npoints kdtree) ;; index of Next Free slot
	with indexvec-left = (kdtree-index-left-vec kdtree)
	with indexvec-right = (kdtree-index-right-vec kdtree)
	with flagvec = (kdtree-flag-vec kdtree)
	with ir-vec = (kdtree-ir-vec kdtree)
	with r-vec = (kdtree-r-vec kdtree)
	with ndim-1 = (1- ndim)
	for idepth of-type fixnum from 0 
	;; the dimension we're splitting on
	for idim = 0 then (if (= idim ndim-1) 0 (1+ idim)) ;; faster than mod
	  ;;initially   (incf (kdtree-npoints kdtree) 1)
	do
	   ;(incf *nops*)
	   (cond
	     ;; this node is empty, so put our data here
	     ((= (aref flagvec k) +empty+)
	      (%insert-into-kd v object kdtree k k idepth idim)
	      (return k))
	     ;;
	     ;; splitting in left direction - note indirection with ir-vec
	     ((<= (aref v idim) (aref r-vec (aref ir-vec k) idim))
	      (cond
		;; no next left element, so insert it at k+1
		((= (aref indexvec-left k) +end+)
		 (setf (aref indexvec-left k) knf) ;; next free element
		 (%insert-into-kd v object kdtree knf knf idepth
				  (if (= idim ndim-1) 0 (1+ idim)))
		 (return knf))
		;; else traverse further
		(t
		 (setf k (aref indexvec-left k)))))
	     ;;
	     ;; splitting in right direction
	     (t
	       (cond
		;; no next left element, so insert it at k+1
		((= (aref indexvec-right k) +end+)
		 (setf (aref indexvec-right k) knf) ;; next free element
		 (%insert-into-kd v object kdtree  knf knf idepth
				  (if (= idim ndim-1) 0 (1+ idim)))
		 (return knf))
		;; else traverse further
 		(t
		 (setf k (aref indexvec-right k))))))))


(defun insert-2d (kdtree x y object &key (vec nil) (defer nil))
  "Insert X,Y into 2d KDTREE."
  (declare (type kdtree kdtree)
	   (type real x y)
	   (type (or null kdfltvec) vec))
  (let ((vec (or vec (make-kdfltvec 2))))
    (setf (aref vec 0) (coerce x 'kd-float))
    (setf (aref vec 1) (coerce y 'kd-float))
    (insert-vector kdtree vec object :defer defer)))

(defun insert-3d (kdtree x y z object &key (vec nil) (defer nil))
  "Insert X,Y,Z into 3d KDTREE."
  (declare (type kdtree kdtree)
	   (type real x y z)
	   (type (or null kdfltvec) vec))
  (let ((vec (or vec (make-kdfltvec 3))))
    (setf (aref vec 0) (coerce x 'kd-float))
    (setf (aref vec 1) (coerce y 'kd-float))
    (setf (aref vec 2) (coerce z 'kd-float))
    (insert-vector kdtree vec object :defer defer)))

  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defun make-deletion-action (&key (test (constantly t)))
  "Function to generate a function to pass as an :ACTION, to delete
objects according to :TEST.

TEST is an optional function taking the OBJECT inside KDTREE-OBJ-VEC,
returning T if this object is to be deleted.

For example, 
 (MAKE-DELETION-ACTION (LAMBDA (OBJ) (EQ OBJ *MY-OBJ*))) 
will return an action to delete objects EQ to *MY-OBJ*
that are picked up by the search.
and 
 (MAKE-DELETION-ACTION) will return an action to
delete all objects picked up by the KDTREE search."
  (lambda (kdtree inode)
    (declare (type kdtree kdtree)
	     (type index inode)
	     (optimize speed))
    (when (funcall test (aref (kdtree-obj-vec kdtree) (aref (kdtree-ir-vec kdtree) inode)))
      (setf (aref (kdtree-obj-vec kdtree) (aref (kdtree-ir-vec kdtree) inode)) 'deleted-object))))

