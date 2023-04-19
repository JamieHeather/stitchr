Stitching immunoglobulins
=========================

BCR and antibodies are of course produced in a similar manner to TCRs, and thus the conceptual framework applies equally to those sequences. However, the existence of somatic hypermutation, far greater structural- and allelic-polymorphism, and more complicated constant region biology in the IG loci makes ``stitchr`` more difficult to apply. While we are not officially supporting the use of ``stitchr`` to make immunoglobulin sequences, we have made it possible for users to do so if they wish.

In order to achieve this, users can try using `this modified version of stitchrdl <https://gist.github.com/JamieHeather/bb86ea2c5968bf93a6c1df59c8849435>`_ to download formatted human immunoglobulin sequences from IMGT. This will create an additional species folder called 'HUMAN-IG' in the stitchr data folder, wherever that's installed. Note that we strongly recommend thoroughly checking the sequences downloaded in this manner. Alternatively users may wish to curate their own germline datasets as per the :ref:`data-formatting-label` section, and move it to the stitchr data directory (``stitchr -dd``) themselves.

These can then be used like so:

.. code:: bash

   stitchr -v IGHV3-30-3*01 -j IGHJ4*02 -cdr3 CARLSPAGGFFDYW -c IGHM*01 -s HUMAN-IG -n JQ304252
   stitchr -v IGHV4-61*01 -j IGHJ3*02 -cdr3 CARITGDRGAFDIW -c IGHD*01 -s HUMAN-IG -n AF262208
   stitchr -v IGHV1-69*01 -j IGHJ3*02 -cdr3 CAREVVPTFRENAFDIW -c IGHG1*01 -s HUMAN-IG -n MW177368
   stitchr -v IGHV4-59*01 -j IGHJ5*02 -cdr3 CARGISWFDPW -c IGHE*01 -s HUMAN-IG -n DQ005305
   stitchr -v IGKV3-20*01 -j IGKJ5*01 -cdr3 CQQYGTSRPITF -c IGKC*01 -s HUMAN-IG -n BC032451
   stitchr -v IGLV1-47*01 -j IGLJ3*02 -cdr3 CAAWDDSLSGWVF -c IGLC2*01 -s HUMAN-IG -l IGLV1-47*02 -n AB064224

However for the reasons stated above we recommend using caution when applying ``stitchr`` to these loci: long read sequencing (both into the V and the C) and liberal use of the ‘seamless’ setting is recommended.

Note that the default form of the IGH constant regions supplied when using just the gene+allele is the secreted form: (where available) the membrane bound form is produced by appending ’_M’. E.g. use ``IGHM*01`` for the secreted form, and ``IGHM*01_M`` for the membranous.
