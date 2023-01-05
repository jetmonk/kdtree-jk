

(asdf:defsystem kdtree-jk
  :author "jetmonk@gmail.com"
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
  :license "MIT"
  :depends-on (kdtree-jk)
  :components
  ((:file "kdtree-jk-latlon")))


(asdf:defsystem kdtree-jk/tests
  :author "jetmonk@gmail.com"
  :license "MIT"
  :depends-on (kdtree-jk kdtree-jk/latlon)
  :components
  ((:file "kdtree-jk-tests")
   (:file "kdtree-jk-tests-latlon")))
