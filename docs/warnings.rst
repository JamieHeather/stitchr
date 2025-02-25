
``stitchr`` warnings
~~~~~~~~~~~~~~~~~~~~

As stitchr runs it will often generate a number of warning strings which describe properties of the rearrangement in question. While some of these sound bad, they are often perfectly fine, usually reflecting either a potential minor annotation or data input error, or a rearrangement with some slightly-uncommon recombination parameters.

Some of the warnings most commonly seen are explained below:

* *"Note: while a C-terminal CDR3:J germline match has been found, it was only the string 'F'*" (or something similar)

    - The way that ``stitchr`` figures out how far to trim the germline regions is to effectively take the ends of the CDR3 and ask 'how much of this could be encoded by the V and the J?'. It does that by taking the amino acid sequence of say the 3' end of the CDR3, and looking for that amino acid sequence in the translation of the J gene. If it doesn't find one, it trims off one of the N terminal residues and then looks again.

    * As most J genes end their CDR3s with a conserved phenylalanine, 'F' is often the minimum amount that can be matched.

    * This can present technical problems in cases where say more than one F occurs in a J gene, or if people/tools annotate their CDR3 junctions improperly (e.g. perhaps ending in an F that's not the conserved F, which is an error that occurs more frequently in J genes that have an additional F immediately upstream of the FGXG motif).

    * However even if ``stitchr`` has detected the corrected conserved F residue, this position occurs relatively far into the J gene, past the point where most J genes would get deleted to during recombination - which is why ``stitchr`` flags this up as a potential problem.

    - Usually this indicates the incorrect J allele being provided, so it's looking for the wrong amino acid sequence. E.g. a common example is seen in VDJdb data of TCRs using TRAJ24, as many of the \*01 allele calls are likely actually supposed to be \*02. As these two sites differ by the two residues immediately upstream of the conserved F (below) misannotation prevents a typical length J motif being found:

::

        TTDSWGKFEFGAGTQVVVTP = TRAJ24*01
        TTDSWGKLQFGAGTQVVVTP = TRAJ24*02
        *******—-***********

* *Warning: TR[AB]JXX*YY has a 'low confidence' CDR3-ending motif.*

    * This warning is produced when ``stitchr`` detects that a non-classical FGXG motif is used in a particular J gene. This doesn't necessarily mean that this J (or TCRs using it) are non-functional, but as this is atypical it gets flagged for users to apply appropriate caution - especially in species for which there may be less data.

    * E.g. this commonly occurs with TRAJ35 in human data, as it uses a 'CGSG' motif, yet is seemingly functional.


* *Unable to locate N terminus of CDR3 in V gene correctly. Please ensure sequence plausibility.*

    * This is similar to the V gene equivalent of the first J gene problem above, in that it occurs when ``stitchr`` cannot find the expected 2nd-Cys motif that should form the beginning of the CDR3 function.

    * This usually occurs as a failure of whatever VDJ annotation software was used, typically partial or truncated CDR3 sequences.


* *Error: a [LEADER/VARIABLE/JOINING/CONSTANT] sequence region has not been found for gene in the IMGT data for this chain/species. Please check your TCR and species data.*

    * It's a requirement of ``stitchr`` that it has to have an entry for every part of every gene (leader/V/J/constant) that it's being asked to stitch. This error occurs when you've asked ``stitchr`` to use a gene that it doesn't have all of the necessary sequences for.

    * This error most frequently occurs for leader regions (especially in non-human species), even when the corresponding V gene is included. In these circumstances you can still get ``stitchr`` to produce something if you don't care about the 'true' sequence, by using a different V gene or allele's leader sequence, or by adding the missing sequence to the appropriate file in the data directory. Alternatively, if you don't care about expressing the stitched sequence, you can add any arbitrary DNA sequence (e.g. '``ATG``', or even '``nnn``' will work).
