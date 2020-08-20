(in-package :load-a-pdb)

;;; Make this look like tirun-jupyter.lisp


(defmacro def-dyn-widget (name initial-value)
  (let ((var-name (gensym)))
    `(progn
       (defparameter ,var-name nil)
       (defun ,name ()
         (or ,var-name
             (setf ,var-name ,initial-value))))))


(def-dyn-widget button-style (make-instance 'w:button-style :button-color "aquamarine"))


(defclass app (jupyter-widgets:has-traits)
  ((receptor-string
    :accessor receptor-string
    :initform nil
    :trait :string)
   (loaded-ligands
    :accessor loaded-ligands
    :initform nil
    :trait :list)
   (selected-ligands
    :accessor selected-ligands
    :initform nil
    :trait :list)
   (ligand-widget
    :accessor ligand-widget
    :initform nil)
   (all-nodes
    :accessor all-nodes
    :initform nil
    :trait :list)
   (all-edges
    :accessor all-edges
    :initform nil
    :trait :list)
   (smirks
    :accessor smirks
    :initform ""
    :trait :string)
   (distributor
    :accessor distributor
    :initform "s103.thirdlaw.tech"
    :trait :string)
   (job-name
    :accessor job-name
    :initform "default"
    :trait :string)
   (submit-stream
    :accessor submit-stream
    :initform nil
    :trait :stream)
   (cyto-observe
    :accessor cyto-observe
    :initform nil
    :trait :list)
   (cyto-widget
    :accessor cyto-widget
    :initform nil
    :trait :list)
   (smirk-pattern
    :accessor smirk-pattern
    :initform nil
    :trait :string)



   )  
  (:metaclass jupyter-widgets:trait-metaclass))

