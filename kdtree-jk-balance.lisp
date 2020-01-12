

(in-package kdtree-jk)



(defun balance-kdtree (kdtree)
  "Balances a KDTREE by ensuring that each node splits the remaining set in 2.
Also removes any deleted objects (objects with object-vec set to 'DELETED-OBJECT)"
  (declare (type kdtree kdtree)
	   (optimize speed))

  ;; if any deleted objects (obj-vec[j]='deleted-object) then
  ;; compact the arrays to remove them
  (%remove-deleted-objects kdtree)
  
  (let* ((np (kdtree-npoints kdtree))
	 (ndim (kdtree-ndim kdtree))
	 (flag-vec (kdtree-flag-vec kdtree))
	 (index-left-vec (kdtree-index-left-vec kdtree))
	 (index-right-vec (kdtree-index-right-vec kdtree))
	 (ir-vec (kdtree-ir-vec kdtree)) ;; this gets re-arranged, and replaced by ir-vec2
	 (ir-vec2 (make-indexvec (kdtree-npoints-total kdtree))) ;; the new vector constructed below
	 (sd-vec (kdtree-sd-vec kdtree)))

    (declare (type indexvec ir-vec2))
     
    (loop for j below (kdtree-npoints-total kdtree)
	  do (setf (aref flag-vec j) +empty+)
	     (setf (aref index-left-vec j) +end+)
	     (setf (aref index-right-vec j) +end+)
	     (setf (aref sd-vec j) 0))
    (setf (kdtree-idepth kdtree) 0)
    (setf (kdtree-avg-depth kdtree) 0d0)
    (setf (kdtree-npoints kdtree) 0)

    ;; in the recursion, (kdtree-npoints kdtree) is the index of our current position
    (labels ((recurse-balance (j1 j2 idim idepth)
	       (declare (type index j1 j2 idepth)
			(type dimnum idim))
	       (let* ((j (if (= j1 j2)
			     j1 ;; last node
			     (partition-ir-vec-around-median kdtree j1 j2 idim)))
		      (idim-next (if (= idim (1- ndim)) 0  (1+ idim)))
		      (k (kdtree-npoints kdtree))) ;; k is the index of this node 
		 (declare (type index j k)
			  (type dimnum idim-next))

		 (incf (kdtree-npoints kdtree))
		 (setf (aref (kdtree-flag-vec kdtree) k) +filled+)
		 (setf (aref (kdtree-sd-vec kdtree) k) idim)
		 (setf (aref ir-vec2 k) (aref ir-vec j))
		 (setf (kdtree-idepth kdtree) (max idepth (kdtree-idepth kdtree)))
		 (let ((f (/ 1d0 (kdtree-npoints kdtree)))) ;;fraction of points this point is
		   (setf (kdtree-avg-depth kdtree)
			 (+ (* (- 1d0 f) (kdtree-avg-depth kdtree))
			    (* f idepth))))
		 ;; fill left tree if there are points remaining
		 (when (<= j1 (- j 1))
		   (setf (aref index-left-vec k)
			 (recurse-balance j1 (- j 1)
					  idim-next
					  (1+ idepth))))
		 ;; fill right tree if there are points remaining
		 (when (>= j2 (+ j 1))
		   (setf (aref index-right-vec k)
			 (recurse-balance (+ j 1) j2
					  idim-next
					  (1+ idepth))))
		 ;; return the current index of this node to put into index-XXXX-vec one level higher
		 k)))

      (recurse-balance 0 (1- np) 0 0)
      (setf (kdtree-ir-vec kdtree) ir-vec2)
      (setf (kdtree-needs-balancing kdtree) nil)
      t)))
     
	      
