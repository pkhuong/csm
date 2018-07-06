(defun scale-for (width)
  (check-type width (real (0) 1))
  (loop with scale = 1d0 while (<= (* scale width) 1d0)
     do (setf scale (* scale 10d0))
     finally (return scale)))

(defun range-and-scale (min max &key (subticks 5))
  (check-type min (real 0 1))
  (check-type max (real 0 1))
  (when (< max min)
    (rotatef min max))
  (let* ((width (max (- max min) 1d-18))
         (scale (scale-for width))
         (scaled-min (floor (- (* min scale) 0.0d0)))
         (scaled-max (ceiling (+ (* max scale) 0.0d0))))
    (values (* subticks scale)
            (* subticks (max scaled-min 0))
            (* subticks (min scaled-max scale)))))

(defun denormalized-beta-pdf (a b x)
  (when (< (min x (- 1d0 x)) double-float-epsilon)
    (return-from denormalized-beta-pdf double-float-negative-infinity))
  (let ((log-x (log (float x 1d0)))
        (log-1mx (sb-kernel:%log1p (- x))))
    (+ (* log-x (1- a)) (* log-1mx (1- b)))))

(defun beta-pdf-points-1 (a b scale min max max-value)
  (flet ((point (x)
           (list x (denormalized-beta-pdf a b (/ x scale))))
         (high-resolution (x)
           (let* ((points (list x (+ x 0.5d0) (1+ x)))
                  (pdfs (mapcar (lambda (x)
                                  (denormalized-beta-pdf a b (/ x scale)))
                                points)))
             (ceiling (- (exp (- (reduce #'max pdfs) max-value))
                         (exp (- (reduce #'min pdfs) max-value)))
                      0.01))))
    (loop for i from min to max
       collect (point i)
       when (< i max)
       append (loop with subscale = (min (high-resolution i) 10)
                 for j from 1 below subscale
                 collect (point (+ i (/ j (float subscale 1d0))))))))

(defun beta-pdf-mode-density-estimate (a b min max)
  (let ((mode (if (or (< a 1) (< b 1))
                  (* .5d0 (+ min max))
                  (/ (1- a) (+ a b -2d0)))))
    (reduce #'max (list* min max (if (< min mode max)
                                     (list mode)
                                     '()))
            :key (lambda (x)
                   (denormalized-beta-pdf a b x)))))

(defun beta-pdf-points (a b scale min max)
  (beta-pdf-points-1 a b scale min max
                     (beta-pdf-mode-density-estimate a b
                                                     (/ min scale)
                                                     (/ max scale))))

(defun rescaled-pdf-points (points)
  (let ((max (reduce #'max points :key #'second
                     :initial-value double-float-negative-infinity)))
    (mapcar (lambda (pair)
              (destructuring-bind (x v) pair
                (list x (exp (if (eql v double-float-negative-infinity)
                                 v
                                 (- v max))))))
            points)))

(defun closed-path-string (points)
  (let ((*read-default-float-format* 'double-float))
    (unless (null points)
      (format nil "M ~A 0 ~{L ~{~,4F ~,4F~} ~} L ~A 0 Z"
              (first (first points))
              points
              (first (first (last points)))))))

(defun call-with-svg-viewbox (out x-min x-max y-min y-max function)
  (let ((*read-default-float-format* 'double-float))
    (format out "<svg viewBox=\"~A ~A ~A ~A\" ~
                 xmlns=\"http://www.w3.org/2000/svg\">~%"
            x-min y-min (- x-max x-min) (- y-max y-min))
    (funcall function out)
    (format out "</svg>")))

(defun call-with-group-class (out class function)
  (format out "<g class='~A'>~%" class)
  (funcall function out)
  (format out "</g>~%"))

(defun draw-distribution (out points max-height)
  (format out "<path d=\"~A\"
               fill=\"red\" fill-opacity=\"0.25\" ~
               stroke=\"transparent\" stroke-width=\"0.01\" ~
               transform=\"scale(~,4F, ~,4F) translate(~,4F, 0)\" />~%"
          (closed-path-string points)
          (/ 1d0
             (- (first (first (last points))) (first (first points))))
          (- (float max-height 1d0))
          (- (first (first points)))))

(defun draw-ticks (out x-min x-max subticks)
  (let ((path (with-output-to-string (s)
                (loop for i from 0 upto (- x-max x-min)
                   do (format s "M ~,4F 0 V -2 "
                              (/ i (- x-max x-min)))))))
    (format out "<path d='~A' stroke='white' stroke-width='0.0075' />~%"
            path))
  (let ((path (with-output-to-string (s)
                (format s "M ~A 0 H ~A " x-min x-max)
                (loop for i from x-min to x-max by subticks
                   for counter upfrom 0
                   do (format s "M ~A 0 V 1 " i)))))
    (format out "<path d='~A' stroke='black' stroke-width='0.05' ~
                  transform='scale(~,4F, 0.125) translate(~,4F, 0)'/>~%"
            path
            (/ 1d0 (- x-max x-min))
            (- x-min))))

(defun emit-styles (out styles)
  (format out "<style type='text/css'>~%")
  (loop for (name . size) in styles do
       (format out ".~A { font-size: ~A; }~%" name size))
  (format out "</style>~%"))

(defun write-limits (out scale x-min x-max)
  (format out
          "<text x='0' y='0' dx='-0.1' dy='0.2' class='limits'>~,2G</text>~%"
          (/ x-min scale))
  (format out
          "<text x='1' y='0' dx='-1em' dy='0.2' class='limits'>~,2G</text>~%"
          (/ x-max scale)))

(defun draw-bound (out scale x-min x-max bound color height)
  (format out "<path d='M ~,4F 0 V ~,4F' stroke='~A' stroke-width='0.02' />~%"
          (/ (- (* scale bound) x-min) (- x-max x-min))
          (* -1d0 height)
          color))

(defun draw-ci (out scale x-min x-max ci-min ci-max ci-center)
  (let* ((width (- x-max x-min))
         (path-string (with-output-to-string (s)
                        (format s "M ~,4F -0.1 V 0.1 V 0 ~
                                   H ~,4F V -0.05 V 0.05 V 0 ~
                                   H ~,4F V -0.1 V 0.1"
                                (/ (- (* scale ci-min) x-min) width)
                                (/ (- (* scale ci-center) x-min) width)
                                (/ (- (* scale ci-max) x-min) width)))))
    (format out "<path d='~A' stroke='grey' stroke-width='0.02' ~
                  fill='transparent' transform='translate(0, -0.25)' />~%"
            path-string)
    (format out
            "<text x='~,4F' y='-0.31' dx='-1.25em' class='center'>~
             ~,3G</text>~%"
            (/ (- (* scale ci-center) x-min) width)
            ci-center)))

(defun write-symbol (out scale x-min x-max symbol ci-center)
  (let ((width (- x-max x-min)))
    (format out
            "<text x='~,4F' y='-0.075' dx='-0.3em' class='symbol'>~
             ~A</text>~%"
            (/ (- (* scale ci-center) x-min) width)
            symbol)))

(defun output-svg (out bound-1 bound-2 ci-min ci-max ci-center
                   a b &key (subticks 2) symbol)
  (let ((min (min bound-1 bound-2 ci-min ci-max))
        (max (max bound-1 bound-2 ci-min ci-max))
        (ci-center (or ci-center (/ a (+ a b))))
        (symbol (cond (symbol symbol)
                      ((> ci-min (max bound-1 bound-2))
                       "&#x226B;") ; >>
                      ((< ci-max (min bound-1 bound-2))
                       "&#x226A;") ; <<
                      ((and (> ci-min (min bound-1 bound-2))
                            (< ci-max (max bound-1 bound-2)))
                       "&#x2248;") ; ~=
                      ((> ci-min (min bound-1 bound-2))
                       "&#x2265;") ; >=
                      ((< ci-max (max bound-1 bound-2))
                       "&#x2264;") ; <=
                      (t
                       "&#xbf;")))) ; inverted question mark
    (multiple-value-bind (scale x-min x-max)
        (range-and-scale min max :subticks subticks)
      (let ((pdf-points (rescaled-pdf-points
                         (beta-pdf-points a b scale x-min x-max))))
        (call-with-svg-viewbox out -0.2d0 1.2d0 -0.5d0 0.21d0
                               (lambda (out)
                                 (draw-distribution out pdf-points 0.49d0)
                                 (format out "~%")
                                 (draw-ticks out x-min x-max subticks)
                                 (emit-styles out 
                                              '(("base" . "1.5%")
                                                ("limits" . "50%")
                                                ("center" . "75%")))
                                 (call-with-group-class
                                  out "base"
                                  (lambda (out)
                                    (write-limits out scale x-min x-max)
                                    (when (/= bound-1 bound-2)
                                      (draw-bound out scale x-min x-max
                                                  (min bound-1 bound-2)
                                                  "blue" 0.2d0))
                                    (draw-bound out scale x-min x-max
                                                (max bound-1 bound-2)
                                                "red" 0.2d0)
                                    (draw-ci out scale x-min x-max
                                             ci-min ci-max ci-center)
                                    (write-symbol out scale x-min x-max
                                                  symbol ci-center)))))))))

(defun svg-data-url (bound-1 bound-2 ci-min ci-max ci-center
                     a b &key (subticks 5) symbol)
  (format t "data:image/svg+xml;base64,~A"
          (cl-base64:string-to-base64-string
           (with-output-to-string (s)
             (output-svg s bound-1 bound-2 ci-min ci-max ci-center
                         a b :subticks subticks :symbol symbol)))))
