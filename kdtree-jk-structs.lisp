
(in-package kdtree-jk)

(deftype kd-float ()
  'double-float)
(deftype kd-flag ()
  'bit)
(deftype index ()
  '(signed-byte 32))
(deftype dimnum () ;;  integer representing a dimension
  '(unsigned-byte 16))

(defconstant +most-positive-float+
  (if (typep 1d0 'kd-float)
      most-positive-double-float
      most-positive-single-float))

;; an end node indicator (no child to that side)
(defconstant +end+ -1)

;; bit flag in flagvec, if this node has a value yet
(defconstant +filled+ 1)
(defconstant +empty+ 0)

(deftype kdfltvec ()
  '(simple-array kd-float (*)))
(deftype kdfltarr ()
  '(simple-array kd-float (* *)))
(deftype indexvec ()
  '(simple-array index (*))) ;; signed so that -1 can be a flag
(deftype flagvec ()
  '(simple-array kd-flag (*)))
(deftype objvec () ;; array storing the binned object
  `(simple-array t (*)))
(deftype dimvec() ;; array storing 
  `(simple-array dimnum (*)))

(defun make-kdfltvec (n)
  (make-array n :element-type 'kd-float))
(defun make-kdfltarr (nrow ncol)
  (make-array (list nrow ncol) :element-type 'kd-float))
(defun make-indexvec (n &optional (initial-element +end+))
  (make-array n :element-type 'index :initial-element initial-element))
(defun make-flagvec (n &optional (initial-element +empty+))
  (make-array n :element-type 'kd-flag :initial-element initial-element))
(defun make-objvec (n &optional (initial-element nil))
  (make-array n :element-type t :initial-element initial-element))
(defun make-dimvec (n &optional (initial-element 0))
  (make-array n :element-type 'dimnum :initial-element initial-element))

;; bounding box of kdtree, and submitted to search
(defstruct bbox
  (ndim 0    :type (unsigned-byte 32))
  (rmin-vec  (make-kdfltvec 0) :type kdfltvec)
  (rmax-vec  (make-kdfltvec 0) :type kdfltvec))

#|
A node 'i' in the tree has data in r-vec, obj-vec as:

  k=ir-vec[i]
  r-vec[k][1..ndim], obj-vec[k]  ;; actual x,y,z.. data and the stored object
  flag-vec[i]  ;; 1/0 if this node is filled/empty

and the left, right nodes are indexed as

  J=index-xxx-vec[i]    ;; next right,left node
  K=ir-vec[J]           ;; index of redirection
  r-vec[K][1..ndim], obj-vec[K]  ;; actual x,y,z.. data and the stored object


the vectors index-left-vec, index-right-vec have the ir-vec indices
to the left (<=) and right (>) trees.  Ie, if j=index-left-vec[i] 
then the left tree is at j

|#
(defstruct kdtree
  ;; dimensions in tree
  (ndim              0 :type (unsigned-byte 16))
  ;; how many points inserted
  (npoints           0 :type (unsigned-byte 32))
  ;; does it need a balancing? This happens on a /defer insertion
  (needs-balancing   nil :type (or null t))
  ;; how deep does this tree go? should not be >> (log npoints 2)
  (idepth            0 :type (unsigned-byte 20))
  ;; average depth of a point - should be about (log npoints 2)
  (avg-depth         0d0 :type double-float)
  ;; how many points available to use (length of vectors)
  (npoints-total     0 :type (unsigned-byte 32))
  ;;
  ;; 1/0 flag if this node is filled; indexing is same as ir-vec
  (flag-vec  (make-flagvec 0 +empty+) :type flagvec)
  ;; 
  ;; Index of Redirection [IR] from current index-XXXX-vec[i] to
  ;;  rvec[k][idim] as j=index-xxx-vec[i] --> k=irvec[j] -->
  ;;  r-vec[k][1..ndim]; this is to allow subsequent re-sorting or
  ;;  balancing of tree without having to move all of the data in
  ;;  r-vec and obj-vec
  (ir-vec (make-indexvec 0 +end+) :type indexvec)
  ;; indices of nodes to left and right, or -1 for END; indexing is
  ;; same as ir-vec
  (index-left-vec  (make-indexvec 0 +end+)  :type indexvec)
  (index-right-vec (make-indexvec 0 +end+)  :type indexvec)
  ;; type T array for storing the lisp object at this node
  ;; nth object is obj-vec[ir-vec[n]]
  (obj-vec (make-array 0)   :type objvec)
  ;; the dimension number this node is splitting on
  (sd-vec (make-dimvec 0) :type dimvec) ;; mainly for diagnostics 
  ;; rvec[ir-vec[n]][m] is mth dimension of nth point in the tree
  (r-vec  (make-kdfltarr 0 0) :type kdfltarr)
  ;; bounding box
  (bbox (make-bbox) :type bbox))


;; points returned from a KD search
(defstruct kdresult
  (n 0 :type (unsigned-byte 32)) ;; how many results
  (obj-vec  (make-objvec 0) :type objvec) ;; vector of objects
  (index-vec (make-indexvec 0) :type indexvec)
  (dist-vec (make-kdfltvec 0) :type kdfltvec))
	     
		   
(defun build-kdtree (ndim &key (npoints 10))
  "Build a KDTREE with NDIM dimensions and NPOINTS of storage,
initially."
  (make-kdtree
   :ndim ndim
   :npoints-total npoints
   :ir-vec (make-indexvec npoints)
   :index-left-vec (make-indexvec npoints)
   :index-right-vec (make-indexvec npoints)
   :flag-vec (make-flagvec npoints)
   :obj-vec (make-objvec npoints)
   :sd-vec (make-dimvec npoints)
   :r-vec (make-kdfltarr npoints ndim)
   :bbox  (make-bbox :ndim ndim
		     :rmin-vec (make-kdfltvec ndim)
		     :rmax-vec (make-kdfltvec ndim))))
	      
    

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defun describe-node (kdtree k &key (stream t))
  "Describe a node in the KDTREE, to STREAM."
  (cond ((= (aref (kdtree-flag-vec kdtree) k) +empty+)
	 (format stream "EMPTY NODE~%"))
	(t
	 (let ((j (aref (kdtree-ir-vec kdtree) k))
	       (rvec (kdtree-r-vec kdtree))
	       (ileftvec (kdtree-index-left-vec kdtree))
	       (irightvec (kdtree-index-right-vec kdtree))
	       (dsvec (kdtree-sd-vec kdtree))
	       (ndim (kdtree-ndim kdtree)))
	   (flet ((print-node (kk)
		    (when (= kk +end+)
		      (return-from print-node))
		    (let ((jj (aref (kdtree-ir-vec kdtree) kk)))
		      (format stream "  VALUES: [")
		      (loop for idim below ndim
			    do (format stream "~,3F~A" (aref rvec jj idim) (if (= idim (1- ndim)) "" ", ")))
		      (format stream "]~%"))))
	     (format stream "NODE ~A indirected to r-vec[~D][..] and obj-vec[~D]~%" k j j)
	     (format stream "Splitting Dim: ~A~%" (aref dsvec k))
	     (print-node k)
	     (format stream "Left Node: ~A~%" (aref ileftvec k))
	     (print-node (aref ileftvec k))
	     (format stream "Right Node: ~A~%" (aref irightvec k))
	     (print-node (aref irightvec k)))))))
	     
	   
