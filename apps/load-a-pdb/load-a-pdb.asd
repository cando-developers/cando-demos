;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :asdf-user)

(defsystem "load-a-pdb"
  :description "load-a-pdb app"
  :version "0.0.1"
  :author ""
  :licence "Private"
  :depends-on (:cando-jupyter)
  :serial t
;;;  :build-operation asdf:monolithic-compile-bundle-op
;;;  :build-pathname #P"/tmp/tirun.fasb"
  :components
  ((:file "packages")
   (:file "load-a-pdb")
   ))

