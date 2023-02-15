

(asdf:defsystem kdtree-jk
  :author "jetmonk@gmail.com"
  :description "KD-TREE package for searching for nearest neighbors in N points in in M-dimensions in N log(N) time."
  :license "MIT"
  :depends-on ()
  :components
  ((:file "kdtree-jk-package")
   (:file "kdtree-jk-structs"
    :depends-on ("kdtree-jk-package"))
   (:file "kdtree-jk-utils"
    :depends-on ("kdtree-jk-package" "kdtree-jk-structs"))
   (:file "kdtree-jk-balance"
    :depends-on ("kdtree-jk-structs" "kdtree-jk-utils"))
   (:file "kdtree-jk-insert"
    :depends-on ("kdtree-jk-structs" "kdtree-jk-balance"))
   (:file "kdtree-jk-search"
    :depends-on ("kdtree-jk-structs" "kdtree-jk-balance"))))


(asdf:defsystem kdtree-jk/latlon
  :author "jetmonk@gmail.com"
  :description "Extension for KDTREE-JK to search in latitude and longitude by converting 2D latitude and longitude into 3D x,y,z."
  :license "MIT"
  :depends-on (kdtree-jk)
  :components
  ((:file "kdtree-jk-latlon")))


(asdf:defsystem kdtree-jk/tests
  :author "jetmonk@gmail.com"
  :description "Tests for KDTREE-JK comparing results to brute force search."
  :license "MIT"
  :depends-on (kdtree-jk kdtree-jk/latlon)
  :components
  ((:file "kdtree-jk-tests")
   (:file "kdtree-jk-tests-latlon")))
