

(setf *default-pathname-defaults* #P"/Users/meister/Development/cando/extensions/cando/src/tests/dynamics/")
(let* ((top #P"complex.parm7")
       (crd #P"complex.rst7")
       (start (amber:start :topology-file top :coordinate-file crd))
       (min (amber:minimize start))
       (heat (amber:heat min))
       (press (amber:pressurize heat))
       (dynamics (amber:dynamics press :nstlim 100000 :ntwx 1000)))
  (amber:generate-all-code (list start)
                           (list (amber:mdcrd dynamics))
                           :pathname-defaults #P"jobs/"))
