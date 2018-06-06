;;;
;;; Initialize everything
;;;

(progn
  (defparameter *working-directory* #P"~/Development/channel-designs/transmembrane04/")
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
	   (p (geom:add p0 (list (geom:v* dir x)))))
      p))

  (defclass unit-cell ()
    ((unit-s :initarg :unit-s :accessor unit-s)
     (unit-t :initarg :unit-t :accessor unit-t)))

  (defun radians (angle-degrees)
    (* 0.0174533 angle-degrees))
  
  (defun p3-unit-cell (dist &optional (angle-degrees 0.0) (s-dir (geom:vec 1.0 0.0 0.0)))
    (let* ((rot-cell (geom:make-m4-rotate-z (radians angle-degrees)))
	   (s-scaled (geom:v* (geom:vnormalized s-dir) dist))
	   (s-cell (geom:m*v rot-cell s-scaled))
	   (trans (geom:make-m4-rotate-z (* -60.0 0.0174533)))
	   (t-cell (geom:m*v trans s-cell)))
      (make-instance 'unit-cell :unit-s s-cell :unit-t t-cell)))

  (defun p1-unit-cell (s-cell t-cell)
      (make-instance 'unit-cell :unit-s s-cell :unit-t t-cell))

  (defun add-unit-cell (accumulate-agg original-agg unit-cell si ti)
    (let* ((vec-s (geom:v* (unit-s unit-cell) si))
	   (vec-t (geom:v* (unit-t unit-cell) ti))
	   (vec (geom:v+ vec-s vec-t))
	   (trans (geom:make-m4-translate vec))
	   (new-agg (chem:matter-copy original-agg)))
      (chem:apply-transform-to-atoms new-agg trans)
      (chem:map-molecules nil (lambda (m) (chem:add-molecule accumulate-agg m))
			  new-agg)))

  (defparameter *positions* '(
			      ( 0 0 )
			      ( 1 0 )
			      ( 0 1 )
			      ( 1 1 )
			      ( 2 0 )
			      ( 2 1 )
			      ( -1 2 )
			      ( 1 2 )
			      ( 0 2 )
			      ( -1 3 )
			      ( 0 3 )
			      ( 1 3 )
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

(defparameter *working-directory* #P"~/Development/channel-designs/transmembrane04/")

;;;
;;; Load the chemdraw file containing a molecule
;;;
(progn
  (setq *default-pathname-defaults* (pathname *working-directory*))
  (defparameter *cd*
    (with-open-file
        (fin (probe-file "trimemwide.cdxml") :direction :input)
      (chem:make-chem-draw fin)))
  (defparameter *agg* (chem:as-aggregate *cd*)))

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
(energy:minimize *agg*)
(energy:minimize *agg* :restraints-on nil)

(cando:chimera *agg*)

(cando:save-mol2 *agg* "trimemwide.mol2")
(defparameter *agg* (cando:load-mol2 "trimemwide.mol2"))

(defparameter *cell-dimension* 16.5)
(defparameter *array* (build-cell *agg* (p3-unit-cell *cell-dimension*)))

(cando:chimera *array*)

(cando:save-mol2 *array* "trimemwide.mol2")

(energy:minimize *array*)
(cando:chimera *array* *axes*)

(defparameter *array* (cando:load-mol2 "trimemwide.mol2"))

;;; Rename the atoms so that they fit within 4 characters and every
;;; atom has a unique name
(antique:rename-atoms-with-unique-short-names *array*)

;;; Save the aggregate with the short names
(cando:save-mol2 *array* "trimemwide-short.mol2")



(setq *default-pathname-defaults* #P"~/Development/channel-designs/transmembrane04/")

;;; I'm seeing problems in the mol2 files that antechamber generates
;;; so I'm going to have it write out an ac format file and
;;; then save that back as a mol2
;;;
;;; Assign partial charges using antechamber
(progn
  (defparameter *ac-in-pathname*
    (namestring (merge-pathnames #P"trimemwide-short.mol2")))
  (defparameter *ac-out-pathname*
    (namestring (merge-pathnames #P"trimemwide-gas.ac"))))
(progn
  (let ((cmd (format nil "antechamber -i ~a -fi mol2 -o ~a -fo ac -c gas"
		     *ac-in-pathname*
		     *ac-out-pathname*)))
    (ext:system cmd))
  (print "Done antechamber"))





;;; Load the aggregate generated by antechamber
(defparameter *ac-agg* (antechamber:load-antechamber *ac-out-pathname*))

(progn
  (defparameter *carbonyl* (make-cxx-object 'chem:chem-info))
  (chem:compile-smarts *carbonyl* "C1(~O2)-N")
  (defun find-carbonyls (agg)
    (let ((num 0)
	  (max-c 0.0)
	  (min-o 0.0))
      (chem:map-atoms
       nil
       (lambda (a)
	 (incf num)
	 (let ((matchp (chem:matches *carbonyl* a)))
	   (when matchp
	     (let* ((match (chem:get-match *carbonyl*))
		    (atomc (chem:get-atom-with-tag match :1))
		    (atomo (chem:get-atom-with-tag match :2))
		    (chc (chem:get-charge atomc))
		    (cho (chem:get-charge atomo)))
	       (format t "~a(~a) ~a(~a) ~a~%" atomc chc atomo cho (+ chc cho))))))
       agg)
      (format t "Searched ~a atoms~%" num)))
  (defun fix-carbonyls (agg)
    (chem:map-atoms
     nil
     (lambda (a)
       (let ((matchp (chem:matches *carbonyl* a)))
	 (when matchp
	   (let* ((match (chem:get-match *carbonyl*))
		  (atomc (chem:get-atom-with-tag match :1))
		  (atomo (chem:get-atom-with-tag match :2))
		  (chc (chem:get-charge atomc))
		  (cho (chem:get-charge atomo))
		  (new-chc 0.51)
		  (new-cho (+ (- 0.51 ) (+ chc cho))))
	     (chem:set-charge atomc new-chc)
	     (chem:set-charge atomo new-cho)))))
     agg))
  (defun total-charge (agg)
    (let ((total-charge 0.0))
      (chem:map-atoms
       nil
       (lambda (a)
	 (setf total-charge (+ total-charge (chem:get-charge a))))
       agg)
      total-charge)))

(format t "Total charge before fix-carbonyls: ~a~%" (total-charge *ac-agg*))
(fix-carbonyls *ac-agg*)
(format t "Total charge after fix-carbonyls: ~a~%" (total-charge *ac-agg*))

(total-charge *ac-agg*)

;;; Save the aggregate as a mol2 file
(cando:save-mol2 *ac-agg* #P"trimemwide-gas-carbonyl.mol2")

(cando:chimera *ac-agg*)



(progn
  (defparameter *unit-cell* (p3-unit-cell *cell-dimension*))
  (defparameter *x-pbc* (geom:v* (unit-s *unit-cell*) 3))
  (defparameter *y-pbc* (geom:vec 0.0 (* (geom:vy (unit-t *unit-cell*)) 4) 0.0))
  (defparameter *pbc-min* (geom:vec 0.0 0.0 -40.0))
  (defparameter *pbc-max* (geom:vec (geom:vx *x-pbc*) (geom:vy *y-pbc*) 40.0))
  (format t "The pbc limits are: ~a ~a~%" *pbc-min* *pbc-max*)
  )
;; Water box is 50x50x80

(defparameter *water* (cando:load-psf-pdb "waterBox-original.psf" "waterBox-original.pdb"))
(print "done")


(defparameter *water-array* (build-cell *water* (p1-unit-cell (geom:vec 50.0 0.0 0.0) (geom:vec 0.0 50.0 0.0)) '((0 0) (0 1) (1 0) (1 1))))

(cando:chimera *water-array*)

(defclass min-max-box ()
  ((box-min :initarg :box-min :accessor box-min)
   (box-max :initarg :box-max :accessor box-max)))

(defparameter *lower-box* (make-instance 'min-max-box
					 :box-min (geom:vec (geom:vx *pbc-min*)
								(geom:vy *pbc-min*)
								-40)
					 :box-max (geom:vec (geom:vx *pbc-max*)
								(geom:vy *pbc-max*)
								-21)))

(defparameter *upper-box* (make-instance 'min-max-box
					 :box-min (geom:vec 0 0 21)
					 :box-max (geom:vec (geom:vx *pbc-max*)
								(geom:vy *pbc-max*)
								40)))


(defun molecule-in-box (molecule box)
  (let ((center (chem:geometric-center molecule)))
    (and (< (geom:vx (box-min box)) (geom:vx center) (geom:vx (box-max box)))
	 (< (geom:vy (box-min box)) (geom:vy center) (geom:vy (box-max box)))
	 (< (geom:vz (box-min box)) (geom:vz center) (geom:vz (box-max box))))))

(defun isolate-molecules (aggregate box)
  (let ((new-agg (chem:make-aggregate)))
    (chem:map-molecules
     nil
     (lambda (m)
       (when (molecule-in-box m box)
	 (chem:add-molecule new-agg m)))
     aggregate)
    new-agg))


(defparameter *lower-water* (isolate-molecules *water-array* *lower-box*))
(defparameter *upper-water* (isolate-molecules *water-array* *upper-box*))

(cando:chimera *lower-water* *upper-water*)
(cando:chimera *lower-water* *upper-water* *array* *axes*)

(defparameter *water-and-channels* (chem:make-aggregate))

(cando:chimera *agg*)
(cando:chimera *array*)
(cando:chimera *array* *lower-water* *upper-water*)

(cando:chimera *upper-water*)

(defparameter *all0* (cando:merge-into-one-aggregate *ac-agg* *upper-water*))
(defparameter *all* (cando:merge-into-one-aggregate *all0* *lower-water*))

(cando:chimera *all*)

(progn
  (defparameter *types* (chem:get-types *ff*))
  (chem:assign-types *types* *all*))

(cando:save-mol2 *all* #P"trimemwide-water-carbonyl.mol2" :use-sybyl-types nil)

(find-carbonyls *all*)
