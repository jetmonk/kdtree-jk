

(asdf:defsystem kdtree-jk
  :depends-on ()
  :components
  ((:file "kdtree-jk-package")
   (:file "kdtree-jk-structs" :depends-on ("kdtree-jk-package"))
   (:file "kdtree-jk-utils"  :depends-on ("kdtree-jk-package" "kdtree-jk-structs"))
   (:file "kdtree-jk-balance"  :depends-on ("kdtree-jk-structs" "kdtree-jk-utils"))
   (:file "kdtree-jk-insert"  :depends-on ("kdtree-jk-structs" "kdtree-jk-balance"))
   (:file "kdtree-jk-search"  :depends-on ("kdtree-jk-structs" "kdtree-jk-balance"))))

  
(asdf:defsystem kdtree-jk/tests
  :depends-on (kdtree-jk)
  :components
  ((:file "kdtree-jk-tests")))
