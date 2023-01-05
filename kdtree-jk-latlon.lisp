

(defpackage kdtree-jk/latlon
  (:use #:cl #:kdtree-jk)
  (:export
   #:kdtree-latlon #:kdtree-latlon-p
   #:build-kdtree-latlon
   #:insert-latlon
   #:kd-latlon-search-in-angle
   ))

 
(in-package kdtree-jk/latlon)

;; Give the kdtree a different name when used for this purpose
(defstruct (kdtree-latlon (:include kdtree)))

(defun build-kdtree-latlon (&key (npoints 10))
  "Build a 3d kdtree; simple wrapper for kdtree-jk:build-kdtree"
  (build-kdtree 3 :npoints npoints))


;; convert LONGITUDE and LATITUDE angles, in
;; degrees, to x,y,z, and place in VEC3  LATITUDE is zero at equator.
;;
;; Formula is 
;;     X = cos(latitude) sin(longitude)
;;     Y = cos(latitude) cos(longitude)
;;     Z = sin(latitude)
(defun lon-lat-to-xyz-vec (longitude latitude vec3)
  (declare (type kd-float latitude longitude)
	   (type kdtree-jk::kdfltvec vec3) 
           (optimize speed))
  ;; need to optimize safety=0 to rid of all consing
  (locally (declare (optimize (speed 3) (safety 0))
		    (type (real -1e-6 1e6) latitude longitude))
    (let* ((pi/180 #.(/ pi 180))
	   (latrad  (* latitude  pi/180))
           (lonrad  (* longitude pi/180))
           (cos-latitude   (cos  latrad))
           (sin-latitude   (sin  latrad))
           (cos-longitude  (cos  lonrad))
           (sin-longitude  (sin  lonrad)))
      (declare
       (type kd-float latrad lonrad
             cos-latitude sin-latitude cos-longitude sin-longitude))
      (setf (aref vec3 0) (* cos-longitude cos-latitude)
	    (aref vec3 1) (* sin-longitude cos-latitude)
	    (aref vec3 2) sin-latitude)
      vec3)))

(defun insert-latlon (kdtree-latlon longitude latitude object
			     &key vec defer)
  "Insert OBJECT LONGITUDE, LATITUDE into a 3d KDTREE-LATLON.  
VEC is an optional scratch vector of length 3 of type KD-FLOAT.  
DEFER means to to defer insertion into tree until balancing."
  (declare (type kdtree kdtree-latlon)
	   (type kd-float latitude longitude)
	   (type (or null kdtree-jk::kdfltvec) vec)
	   (optimize speed))
  (when vec
    (if (not (= (length vec) 3))
	(error "Provided VEC is not of length 3")))
  (let ((vec (or vec (kdtree-jk::make-kdfltvec 3))))
    (lon-lat-to-xyz-vec longitude latitude vec)
    (insert-vector kdtree-latlon vec object :defer defer)))


(defun kd-latlon-search-in-angle
    (kdtree-latlon longitude latitude angular-distance
     &key vec kdresult sort nearest-point action)
    "Search a 3d KDTREE-LATLON around LONGITUDE,LATIDUDE (degrees) and
return points within ANGULAR-DISTANCE (degrees, along surface of unit
sphere).  VEC is an optional scratch vector of length 3 of type
KD-FLOAT.

As in KDTREE package, KDRESULT is the result structure, sort means to sort
by distance, NEAREST-POINT means to return only single nearest point,
and ACTION is to perform action instead of search.

The distance returned in KDRESULT is modified to represent the along-surface
distance of unit sphere."
    (declare (type kdtree kdtree-latlon)
	     (type kd-float latitude longitude angular-distance)
	     (type (or null kdtree-jk::kdfltvec) vec)
	     (optimize speed))
  
  (when vec
    (if (not (= (length vec) 3))
	(error "Provided VEC is not of length 3")))
  
  (when (not (<= 0 angular-distance 90))
    (error "Angular search distance much be in range [0,90] degrees"))
  
  (let ((vec (or vec (kdtree-jk::make-kdfltvec 3)))
	;; the radius of the sphere (centered on the lon,lat on the
	;; unit sphere) that contains the points
	(radius (* 2 (sin (* #.(/ pi 360) angular-distance)))))
	
    (lon-lat-to-xyz-vec longitude latitude vec)
    (let ((kdres
	    (kd-search-in-radius kdtree-latlon vec  radius
				 :kdresult kdresult :sort sort
				 :nearest-point nearest-point
				 :action action)))
      (declare (type kdresult kdres))
      ;; change the distance from being 3d distance from lon,lat to
      ;; point, to the distance along surface of sphere,
      ;; using dsurf = 2 asin (d/2)
      (loop for i below (kdresult-n kdres)
	    for d of-type kd-float = (aref (kdresult-dist-vec kdres) i)
	    for dsurf of-type kd-float = (* 2 (asin (/ (the (real -1 1) d)
						       2)))
	    do (setf (aref (kdresult-dist-vec kdres) i) dsurf))
      kdres)))
    
  
    
  


