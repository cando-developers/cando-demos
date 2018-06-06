;;;
;;; Initialize everything
;;;

(progn
  (defparameter *working-directory* #P"~/Development/channel-designs/transmembrane03/")
  (require :asdf)
  (let ((central-registry (find-symbol "*CENTRAL-REGISTRY*" :asdf))
        (load-system (find-symbol "LOAD-SYSTEM" :asdf)))
    (setf (symbol-value central-registry)
          (cons (translate-logical-pathname "cando:lisp;cando;")
                (symbol-value central-registry)))
    (funcall load-system "cando"))
  (format t "Done initialization pid = ~a~%" (getpid))
  (setq *print-array* t)
  (setq *default-pathname-defaults*
        #P"~/Development/channel-designs/transmembrane03/")
  (let ((*default-pathname-defaults* (probe-file #P"cando:data;force-field;")))
    (defparameter *parms*
      (let ((parms (chem:make-read-amber-parameters)))
	(with-open-file (fin "ATOMTYPE_GFF.DEF" :direction :input)
	  (chem:read-types parms fin))
	(with-open-file (fin "gaff.dat" :direction :input)
	  (chem:read-parameters parms fin)
	  parms)))
    (defparameter *ff* (chem:get-force-field *parms*)))
  (format t "Done with initialization~%"))


(progn
  (defun line (x p0 p1)
    (let* ((dir (geom:sub p1 p0))
	   (p (geom:add p0 (list (geom:times-scalar dir x)))))
      p))

  (defclass unit-cell ()
    ((unit-s :initarg :unit-s :accessor unit-s)
     (unit-t :initarg :unit-t :accessor unit-t)))

  (defun radians (angle-degrees)
    (* 0.0174533 angle-degrees))
  
  (defun p3-unit-cell (dist &optional (angle-degrees 0.0) (s-dir (geom:make-v3 1.0 0.0 0.0)))
    (let* ((rot-cell (geom:make-m4-rotate-z (radians angle-degrees)))
	   (s-scaled (geom:times-scalar (geom:v3-normalized s-dir) dist))
	   (s-cell (geom:mul-v3 rot-cell s-scaled))
	   (trans (geom:make-m4-rotate-z (* 60.0 0.0174533)))
	   (t-cell (geom:mul-v3 trans s-cell)))
      (make-instance 'unit-cell :unit-s s-cell :unit-t t-cell)))

  (defun add-unit-cell (accumulate-agg original-agg unit-cell si ti)
    (let* ((vec-s (geom:times-scalar (unit-s unit-cell) si))
	   (vec-t (geom:times-scalar (unit-t unit-cell) ti))
	   (vec (geom:+ vec-s vec-t))
	   (trans (geom:make-m4-translate vec))
	   (new-agg (chem:matter-copy original-agg)))
      (chem:apply-transform-to-atoms new-agg trans)
      (chem:map-molecules nil (lambda (m) (chem:add-molecule accumulate-agg m))
			  new-agg)))

  (defparameter *positions* '(
			      ( -1 0)
			      ( 0 -1)
			      ( 0 0 )
			      ( 1 0 )
			      ( 0 1 )
			      ( 1 -1 )
			      ( -1 1 )
			      ))

  (defun build-cell (orig-agg unit-cell &optional (positions *positions*))
    (let ((agg (chem:make-aggregate)))
      (dolist (pos positions)
	(let ((ps (first pos))
	      (pt (second pos)))
	  (add-unit-cell agg orig-agg unit-cell ps pt)))
      agg))

  (defun energy-for-agg (agg)
    (let ((ef (chem:make-energy-function agg *ff*)))
      (chem:calculate-energy ef)))

  (defun quick-minimize (energy-func)
    (let ((minimizer (chem:make-minimizer :energy-function energy-func))
	  (restraint-term (chem:get-anchor-restraint-component energy-func)))
      (chem:disable restraint-term)
      (cando:configure-minimizer minimizer
				 :max-steepest-descent-steps 5
				 :max-conjugate-gradient-steps 5
				 :max-truncated-newton-steps 0)
      (chem:enable-print-intermediate-results minimizer)
      #+(or)(chem:set-option energy-func 'chem:nonbond-term nil)
      #+(or)(cando:minimize-no-fail minimizer)
      (chem:set-option energy-func 'chem:nonbond-term t)
      (cando:minimize-no-fail minimizer)))

  (defun check-interactions (agg)
    (let ((energy-function (chem:make-energy-function agg *ff*)))
      (chem:check-for-beyond-threshold-interactions energy-function)))
  
  (defun sample (agg)
    (let ((energy-function (chem:make-energy-function agg *ff*)))
      (quick-minimize energy-function)
      (values (chem:calculate-energy energy-function) (chem:check-for-beyond-threshold-interactions energy-function))))
  
  (defun search-p3-dist (orig-agg &key dist-start dist-end (dist-inc 2.0)
				    rot-start rot-end (rot-inc (* 20.0 0.0174533)))
    (let (results)
      (loop for dist from dist-start to dist-end by dist-inc
	 do (loop for rot from rot-start to rot-end by rot-inc
	       for unit-cell = (p3-unit-cell dist rot)
	       for agg = (build-cell orig-agg unit-cell *positions*)
	       do (progn
		    (format t "Searching distance: ~a rotation: ~a~%" dist rot)
		    (multiple-value-bind  (energy beyond-threshold)
			(sample agg)
		      (push (list energy beyond-threshold dist rot) results)
		      (format t "ENERGY: ~a beyond-threshold: ~a  @ distance: ~a rotation: ~a~%" energy beyond-threshold dist rot )))))
      results))
  )


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Work with a molecule
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *working-directory* #P"~/Development/channel-designs/transmembrane03/")

;;;
;;; Load the chemdraw file containing a molecule
;;;
(progn
  (setq *default-pathname-defaults* (pathname *working-directory*))
  (defparameter *cd*
    (with-open-file
        (fin (probe-file "trimembrane02.cdxml") :direction :input)
      (chem:make-chem-draw fin)))
  (defparameter *agg* (chem:as-aggregate *cd*))
  (defparameter *axes* (merge-pathnames #P"xyz50.wrl" *working-directory*)))

;;;
;;; Gather the atoms that are stereocenters into a vector
;;; and assign them all "S" stereochemistry (#b111...)
(progn
  (defparameter *stereocenters*
    (sort (cando:gather-stereocenters *agg*) #'string< :key #'chem:get-name))
  (cando:set-stereoisomer-func *stereocenters* (constantly :S) :show t)
  (cando:set-stereoisomer-func (select:atoms-with-property *agg* :lead) (constantly :R) :show t)
  (dolist (a (select:atoms-with-property *agg* :stereo))
    (chem:set-configuration a (getf (chem:get-properties a) :stereo))
    (format t "Set ~a to ~a~%" (chem:get-name a) (chem:get-configuration a))))

(defparameter *fix-atoms*
  (sort (select:atoms-with-property *agg* :fix) #'string<
	:key (lambda (a) (string (getf (chem:get-properties a) :fix)))))

*fix-atoms*
(dolist (a *fix-atoms*)
  (format t "properties: ~a~%" (chem:get-properties a)))

(let ((macrocycle-points (anchor:circle-points 10 3 :z 20))
      (top-trifold-points (anchor:circle-points 6 3  :z 15))
      (bottom-trifold-points (anchor:circle-points 10 6 :z -20)))
  (defparameter *points* (append macrocycle-points top-trifold-points bottom-trifold-points))
;;; anchor the :fix atoms to *fixed-points*
  (anchor:on-points *fix-atoms* *points*)
  )

(energy:setup-amber)
(cando:jostle *agg* 40)
(cando:jostle *agg* 4)
(energy:minimize *agg*)
(energy:minimize *agg* :restraints-on nil)

(cando:chimera *agg*)

(cando:save-mol2 *agg* "trichannel.mol2")


(defparameter *array* (build-cell *agg* (p3-unit-cell 20.0)))

(cando:save-mol2 *array* "trichannel-array.mol2")

(cando:chimera *array*)

(defun load-membrane ()
  "Load the membrane from the psf/pdb files 
and name each of the molecules
within it based on whether they 
are solvent (<= 3 atoms) or lipid."
  (let ((agg-membrane
	 (cando:load-psf-pdb "POPC36.psf" "POPC36.pdb")))
    (chem:map-molecules
     nil
     (lambda (m)
       (if (<= (chem:number-of-atoms m) 3)
	   (chem:set-name m :solvent)
	   (chem:set-name m :lipid)))
     agg-membrane)
    agg-membrane))

(defun overlap-func (molecule)
  "* Arguments
- molecule :: A molecule.
* Description
Return the minimum distance for overlap of other atoms with this molecule.
If the molecule is a solvent molecule return 3.0.
If the molecule is a lipid molecule return 0.6."
  (if (eq (chem:get-name molecule) :solvent)
      3.0
      0.6))


(defparameter *agg-membrane* (load-membrane))

;; Remove the bond between hydrogens of TIP3P waters
;; They confuse leap
(let (hh-bonds)
    (chem:map-bonds nil (lambda (a1 a2 o)
			  (when (and (eq (chem:get-element a1) :H)
				     (eq (chem:get-element a2) :H))
			    (push (list a1 a2) hh-bonds)))
		    *agg-membrane*)
    (dolist (atoms hh-bonds)
      (chem:remove-bond-to (first atoms) (second atoms))))

;; Remove molecules from *agg-membrane* that overlap with atoms in *agg-channel*
(cando:remove-overlaps *agg-membrane* *array* :distance-function #'overlap-func)

(cando:chimera *agg-membrane*)
;; Save the hollowed out membrane for reference
(cando:save-mol2 *agg-membrane* #P"hollowed-out-membrane.mol2" :use-sybyl-types t)
(print "Done")

;; Merge the molecules from both aggregates into a new aggregate
(defparameter *agg* (cando:merge-into-one-aggregate *agg-membrane* *ac-agg*))

;; Take a look at the new aggregate
(cando:chimera *agg*)


;; Assign atom types to aggregate

(progn
  (defparameter *types* (chem:get-types *ff*))
  (chem:assign-types *types* *agg*))

;; Save the new aggregate
(cando:save-mol2 *agg* #P"chanmon5-in-membrane.mol2")














(defparameter *results*
  (search-p3-dist *agg* :dist-start 10 :dist-end 16 :dist-inc 2.0
		:rot-start 20 :rot-end 60 :rot-inc 20))

(defparameter *ef* (chem:make-energy-function *a* *ff*))
(chem:calculate-energy *ef*)


  (defparameter *tiled* (tile *agg* *positions*))

  (cando:chimera *axes* *tiled*)
  )

(cando:chimera *axes* *tiled*)
(energy:minimize *tiled* :restraints-on nil)
(print "Hello")

(progn
  (defparameter *a1* (chem:matter-copy *agg*))
  (chem:apply-transform-to-atoms *a1* *t*)
  (chem:apply-transform-to-atoms *a1* *m*)
  (defparameter *a2* (chem:matter-copy *a1*))
  (chem:apply-transform-to-atoms *a2* *m*)
  (defparameter *a3* (chem:matter-copy *a2*))
  (chem:apply-transform-to-atoms *a3* *m*)
  (defparameter *a4* (chem:matter-copy *a3*))
  (chem:apply-transform-to-atoms *a4* *m*)
  (defparameter *a5* (chem:matter-copy *a4*))
  (chem:apply-transform-to-atoms *a5* *m*)
  (defparameter *a6* (chem:matter-copy *a5*))
  (chem:apply-transform-to-atoms *a6* *m*)
  )

(progn
  (cando:chimera *a1*)
  (cando:chimera *a2*)
  (cando:chimera *a3*)
  (cando:chimera *a4*)
  (cando:chimera *a5*)
  (cando:chimera *a6*))

(line 0)
	








(defparameter *agg* (cando:load-mol2 "trifold3-built.mol2"))

(apropos "sexp")
(cando:as-string *agg*)
(let ((*print-readably* t)
      (*print-pretty* nil)
      (*print-circle* t))
  (core:encode  *agg*))


(defparameter *lead-atoms* (select:atoms-with-property *agg* :lead))
(defparameter *apos* (chem:get-position (first *lead-atoms*)))
(defparameter *bpos* (chem:get-position (second *lead-atoms*)))
(defparameter *cpos* (chem:get-position (third *lead-atoms*)))

(defparameter *to-origin* ((geom:times-scalar *apos* -1))
(defparameter 
(energy:minimize *agg* :restraints-on nil)


(cando:jostle *agg* 1)
(energy:minimize *agg*)
(energy:minimize *agg* :restraints-on nil)

(cando:chimera *agg*)

(

(load "~/Downloads/quicklisp-5.lisp")
(load #P"/Users/meister/quicklisp/setup.lisp")
(print "Done")

(cando:save-mol2 *agg* "trifold3-built.mol2")


(defparameter *a* (cando:load-mol2 "trifold3-built.mol2"))

(cando:chimera *a*)

